"""
    Adjoint Reconstruction 
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc
from firedrake_adjoint import *
from pyadjoint import MinimizationProblem, ROLSolver
from pyadjoint.tape import no_annotations, Tape, set_working_tape
import ROL
#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

#logging.set_log_level(1)
#logging.set_level(1)

# Geometric Constants:
y_max = 1.0
x_max = 1.0

#   how many intervals along x/y directions 
disc_n = 100

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 4, 3 
left_id, right_id = 1, 2

# The mesh
mesh = RectangleMesh(nx=disc_n, ny=disc_n, Lx=y_max, Ly=x_max)

# spatial coordinates
X  = x, y = SpatialCoordinate(mesh)
h     = sqrt(CellVolume(mesh))
y_abs     = sqrt(y**2)
yhat  = as_vector((0,y)) / y_abs

# Global Constants:
steady_state_tolerance = 1e-7
max_num_timesteps      = 5
target_cfl_no          = 2.5
max_timestep           = 1.00

# Stokes related constants:
Ra                     = Constant(1e5)   # Rayleigh Number

# Below are callbacks relating to the adjoint solutions (accessed through solve).
# Not sure what the best place would be to initiate working tape!
tape = get_working_tape()

# Temperature related constants:
delta_t                = Constant(1e-6) # Initial time-step
kappa                  = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson scheme where theta_ts = 0.5:
theta_ts               = 0.5

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args) 

# Stokes Equation Solver Parameters:
stokes_solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_type': 'mumps',
    'mat_type': 'aij'
}

# Temperature Equation Solver Parameters:
temperature_solver_parameters = {
        'snes_type': 'ksponly',
        'ksp_converged_reason': None,
        'ksp_monitor': None,
        'ksp_rtol': 1e-2,
        'ksp_type': 'gmres',
        'pc_type': 'sor',
        'mat_type': 'aij'
}


#########################################################################################################
################################## Geometry and Spatial Discretization: #################################
#########################################################################################################

# Set up function spaces - currently using the P2P1 element pair :
V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W    = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Q    = FunctionSpace(mesh, "CG", 2) # Temperature function space (scalar)

# Set up mixed function space and associated test functions:
Z       = MixedFunctionSpace([V, W])
N, M    = TestFunctions(Z)
Y       = TestFunction(Q)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z    = Function(Z)  # a field over the mixed function space Z.
u, p = split(z)     # can we nicely name mixed function space fields?

# Timestepping - CFL related stuff:
ts_func = Function(Q) # Note that time stepping should be dictated by Temperature related mesh.

def compute_timestep(u):
    # A function to compute the timestep, based upon the CFL criterion
    ts_func.interpolate( h / sqrt(dot(u,u)))
    ts_min = ts_func.dat.data.min()
    ts_min = mesh.comm.allreduce(ts_min, MPI.MIN)
    #return min(ts_min*target_cfl_no,max_timestep)
    return 1e-7 


#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Final state, which will be used as reference for minimization, loaded from a file 
final_state = Function(Q, name='RefTemperature')
final_state_file = DumbCheckpoint("./final_state", mode=FILE_READ)
final_state_file.load(final_state, 'Temperature')
final_state_file.close()


# Initial condition
T_ic   = Function(Q, name="T_IC")
#geotherm_checkpoint = DumbCheckpoint("steady_state", single_file=True, mode=FILE_READ)
#geotherm_checkpoint.load(T_ic, 'Temperature')
#geotherm_checkpoint.close()
# Let's start with the final condition
T_ic.project(final_state)

# Set up temperature field and initialise based upon coordinates:
T_old    = Function(Q, name="OldTemperature")

# Having a single hot blob on 1.5, 0.0
T_old.assign(T_ic)

T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old


# ********************** Setup Equations ************************ 
# viscosiy 
mu        = Constant(1.0) # Constant viscosity

# deviatoric stresses
def tau(u): return  mu * (grad(u)+transpose(grad(u)))

# traction field 
def trac(u,p): return dot(tau(u),n) - p*n

# Finalise equations
F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
F_stokes += - (dot(N,yhat)*Ra*T_theta) * dx 
F_stokes += - div(u)* M * dx

# Setting free-slip BC for top and bottom
u_top = Function(V, name='PlateVelocity')
bcu_top_Dir     = DirichletBC(Z.sub(0), u_top, (top_id)) # This will be used when we impose plate velocities
#bcu_top_Free    = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id)) # this is used for free-slip adjoint
bcu_base        = DirichletBC(Z.sub(0).sub(1), 0.0, (bottom_id)) # bottom boundary is always free-slip
bcu_rightleft   = DirichletBC(Z.sub(0).sub(0), 0.0, (left_id, right_id)) # right and left boundaries have only no-outflow condition

# Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# Prescribed temperature for top and bottom
bct_base = DirichletBC(Q, 1.0, bottom_id)
bct_top  = DirichletBC(Q, 0.0, top_id)

# For some reason this only works here!
u, p    = z.split() # Do this first to extract individual velocity, pressure and lagrange multplier fields:
u.rename('Velocity') 
p.rename('Pressure')

# A simulation time to track how far we are
simu_time = 0.0000

# Setting adjoint and forward callbacks, and control parameter
control = Control(T_ic)


# Setting up the file containing all the info
boundary_velocity_fi = DumbCheckpoint("velocity", mode=FILE_READ)


# Documenting the first forward run
fwd_run_fi = File('./adjoint_tic/first_fwd_run.pvd')

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):

    # updating boundary condition
    boundary_velocity_fi.set_timestep(t=simu_time, idx=timestep)
    boundary_velocity_fi.load(u_top, 'Velocity')

    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    solve(F_stokes==0, z, bcs=[bcu_top_Dir, bcu_base, bcu_rightleft], solver_parameters=stokes_solver_parameters)

    # updating time-step based on velocities
    delta_t.assign(compute_timestep(u)) # Compute adaptive time-step

    # Temperature system:
    solve(F_energy==0, T_new, bcs=[bct_base, bct_top], solver_parameters=stokes_solver_parameters)

    fwd_run_fi.write(T_new, u, p)

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log("Timestep Number: ", timestep, " Timestep: ", str('%10.4e' %simu_time))

boundary_velocity_fi.close()

# Initialise functional
beta = 1e-4
# setting up the functional
misfit     = assemble(0.5*(T_new - final_state)**2 * dx)
reg_term   = assemble(0.5*beta*dot(grad(T_ic), grad(T_ic)) * dx)
functional = misfit #+ reg_term


## used for output of adjoint files
#class DoublyIndexedFile:
#    def __init__(self, name):
#        self.name = name
#        self.index0 = 0 
#        self.f = None
#        self._open()
#
#    def _open(self):
#        name = self.name + '_{}.pvd'.format(self.index0)
#        self.f = File(name)
#
#    def next(self):
#        self.index0 += 1
#        self._open()
#
#    def write(self, function):
#        self.f.write(function)
#
#lambda_t_pvd   = DoublyIndexedFile('Visual/lambda_t')
#lambda_u_pvd   = DoublyIndexedFile('Visual/lambda_u')
#lambda_p_pvd   = DoublyIndexedFile('Visual/lambda_p')

# Below are callbacks allowing us to access various field information (accessed through reducedfunctional).
class OptimisationOutputCallbackPost:
    def __init__(self):
        self.iter_idx = 0
        self.opt_file             = File('Optimisation/opt_file.pvd') 
        #self.opt_gradient_file    = File('Visual/opt_gradient.pvd')
        #self.opt_t_ic_file        = File('Visual/opt_temperature_ic.pvd')
        #self.opt_t_final_file     = File('Visual/opt_temperature_final.pvd')
        #self.opt_u_final_file     = File('Visual/opt_velocity_final.pvd')
        #self.opt_p_final_file     = File('Visual/opt_pressure_final.pvd')
        self.grad_copy            = Function(Q, name="Gradient")
        self.T_ic_copy               = Function(Q, name="InitTemperature")
        self.T_tc_copy               = Function(Q, name="FinTemperature")
        self.z_copy               = Function(Z, name="Stokes")
        self.u_copy               = Function(V, name="Velocity")
        self.p_copy               = Function(W, name="Pressure")

    def __call__(self, cb_functional, dj, controls):
        # output current gradient:
        self.grad_copy.assign(dj)
        # output current control (temperature initial condition)
        self.T_ic_copy.assign(controls)
        # output current final state temperature
        self.T_tc_copy.assign(T_new.block_variable.checkpoint)
        #self.opt_t_final_file.write(self.T_copy)
        # output current final state velocity and pressure
        self.z_copy.assign(z.block_variable.checkpoint)
        self.u_copy.assign(self.z_copy.split()[0])
        self.p_copy.assign(self.z_copy.split()[1])
        #self.opt_u_final_file.write(self.u_copy)
        #self.opt_p_final_file.write(self.p_copy)
        #self.opt_gradient_file.write(self.grad_copy)
        #self.opt_t_ic_file.write(self.T_copy)
        self.opt_file.write(self.T_ic_copy, self.T_tc_copy, self.u_copy, self.p_copy, self.grad_copy)
        func_val = assemble((self.T_tc_copy-final_state)**2 * dx) 
        grad_val = sqrt(assemble((self.grad_copy)**2 * dx))

        log(str('Iteration= %i, functional=' %(self.iter_idx)), func_val, ', grad(dj)=', grad_val)
        log('Derivative calculation', self.iter_idx)
        self.iter_idx += 1

        #global lambda_t_pvd, lambda_u_pvd, lambda_p_pvd
        #lambda_t_pvd.next()
        #lambda_u_pvd.next()
        #lambda_p_pvd.next()

class ForwardCallbackPost:
   def __init__(self):
      self.fwd_idx = 0 
   def __call__(self, controls):
      self.fwd_idx +=1 
      log("Starting fwd calculation:", self.fwd_idx)

# Initiate classes for the callbacks
local_cb_post = OptimisationOutputCallbackPost()
eval_cb_pre = ForwardCallbackPost()





# Defining the object for pyadjoint
reduced_functional = ReducedFunctional(functional, control, eval_cb_pre=eval_cb_pre, derivative_cb_post=local_cb_post)

#my_der_l2 =reduced_functional.derivative({'riesz_representation':'l2'});my_der_l2.rename('derl2') 
#my_der_L2 =reduced_functional.derivative({'riesz_representation':'L2'});my_der_L2.rename('derL2') 
#my_der_H1 =reduced_functional.derivative({'riesz_representation':'H1'});my_der_H1.rename('derH1') 
#derivatives_fi = File('derivatives.pvd')
#derivatives_fi.write(my_der_l2, my_der_L2, my_der_H1)
#import sys;sys.exit()

# Set up bounds, which will later be used to enforce boundary conditions in inversion:
T_lb     = Function(Q, name="LB_Temperature")
T_ub     = Function(Q, name="UB_Temperature")
T_lb.assign(0.0)
T_ub.assign(1.0)


### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))

class InitHessian(ROL.InitBFGS):
    def __init__(self, M, direct=True):
        ROL.InitBFGS.__init__(self, M)

        self.solver_params = {
            "ksp_type": "gmres",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "mat_type": "aij",
        }
        # K matrix
        self.M = dot(grad(Y), grad(TrialFunction(Q))) * dx

        # regular linear solver
        self.solver = LinearSolver(assemble(self.M), solver_parameters=self.solver_params)
        

    @no_annotations
    def applyH0(self, Hv, v):
        self.solver.solve(Hv.dat[0], v.dat[0])
        l = Hv.norm() / v.norm()
        Hv.scale(1 / l)
        self.scaleH0(Hv)

    @no_annotations
    def applyB0(self, Bv, v):
        Bv.dat[0] = assemble(dot(grad(Y), grad(v.dat[0])) * dx)
        l = Bv.norm() / v.norm()
        Bv.scale(1 / l)
        self.scaleB0(Bv)
        log('B_0 applied {} times'.format(self.freq))


class StatusReport(ROL.StatusTest):
    def __init__(self, params, vector):
        ROL.StatusTest.__init__(self, params)
        self.vector = vector
        self.T_copy              = Function(Q, name="Temperature")
        self.opt_t_final_file    = File('TEST_2_imposition/opt_temperature_final.pvd')
        self.opt_t_init_file     = File('TEST_2_imposition/opt_temperature_initial.pvd')

    @no_annotations
    def check(self, status):
        log("final state: {}".format(assemble(0.5*(T_new.block_variable.checkpoint - final_state)**2 * dx)))
        log("regularisation: {}".format(assemble(0.5*beta*dot(grad(self.vector.dat[0]), grad(self.vector.dat[0])) * dx)))

        # Writing out final condition
        self.T_copy.assign(T_new.block_variable.checkpoint)
        self.opt_t_final_file.write(self.T_copy)

        ## Writing out initial condition
        self.T_copy.assign(self.vector.dat[0])
        self.opt_t_init_file.write(self.T_copy)
        return ROL.StatusTest.check(self, status)

# This is the classic way
params = {
        'General': {
            #'Secant': {
            #    'Type': 'Limited-Memory BFGS',
            #    'Maximum Storage': 50,
            #    'Line Search': {'Line-Search Method': {'Type': 'Cubic Interpolation'},},}
                'Print Verbosity':1},
        'Step':     {
           'Type': 'Line Search',
           'Line Search': {
                'Descent Method': {'Type': 'Quasi-Newton Method'},
                'Line-Search Method': {'Type': 'Cubic Interpolation'}, 
                'Function Evaluation Limit': 5,
                'Sufficient Decrease Tolerance': 1e-1,
                'Use Previous Step Length as Initial Guess': True,
                            },
        'Status Test': {
            'Gradient Tolerance': 1e-12,
            'Iteration Limit': 50,
                        }
                }
        }


with stop_annotating():    
    # set up ROL problem
    rol_solver = ROLSolver(minp, params)
    sol = rol_solver.solve()



## Orginial parameters by Steph
#params = {
#        'General':  {
#            'Secant':   {
#                'Type': 'Limited-Memory BFGS',
#                'Maximum Storage': 50
#            },
#            'Print Verbosity': 1
#        },
#        'Step': {
#            'Type': 'Line Search',
#            'Line Search':{
#                    'Descent Method': {'Type': 'Quasi-Newton Method'},
#                    'Line Search Method': 'Cubic Interpolation',
#                    'Function Evaluation Limit': 5,
#                    'Sufficient Decrease Tolerance': 1e-1,
#                    'Use Previous Step Length as Initial Guess': True,
#                    }
#                },
#        'Status Test': {
#            'Gradient Tolerance': 1e-12,
#            'Iteration Limit': 2,
#        }
#    }
#
#
#params = ROL.ParameterList(params, "Parameters")
#rol_solver = ROLSolver(minp, params)
##
##init_hessian=0
##if init_hessian == 1:
##    log('Using user-defined initial hessian')
##    rol_secant = InitHessian(5) # maximum storage
##else:
##    log('Using unit matrix for initial hessian')
#rol_secant = ROL.InitBFGS(5) # maximum storage
##rol_step = ROL.TrustRegionStep(rol_secant, params)
##rol_step = ROL.LineSearchStep(params, rol_secant)
#rol_step = ROL.LineSearchStep(params, rol_secant)
###rol_status = ROL.StatusTest(params)
#rol_status = StatusReport(params, rol_solver.rolvector)
#rol_algorithm = ROL.Algorithm(rol_step, rol_status)
##log("functional before going for minimisation: {}".format(assemble(0.5*(T_new - final_state)**2 * dx)))
#
#with stop_annotating():
#    rol_algorithm.run(
#        rol_solver.rolvector,
#        rol_solver.rolobjective,
#        rol_solver.bounds # only if we have a bounded problem
#            )
#
optimal_ic = rol_solver.problem.reduced_functional.controls.delist(rol_solver.rolvector.dat)

# Save the optimal temperature_ic field 
ckpt_T_ic = DumbCheckpoint("T_ic_optimal",\
        single_file=True, mode=FILE_CREATE,\
                           comm=mesh.comm)
ckpt_T_ic.store(optimal_ic)
ckpt_T_ic.close()

