"""
    A finite difference (Taylor's) test for the adjoint
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc
import firedrake.variational_solver as vs
from firedrake_adjoint import * 
from pyadjoint import MinimizationProblem, ROLSolver
from pyadjoint.tape import no_annotations, Tape, set_working_tape
import ROL


#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

x_max = 2.0; y_max = 1.0;
x_dis = 160; y_dis = 80;

# Mesh and associated physical boundary IDs:
mesh    = PeriodicRectangleMesh(x_dis, y_dis, x_max, y_max, direction='x')
bottom_id, top_id = 1, 2  # with PeriodicRectangleMesh: 1,2 for bottom and top

# Global Constants:
steady_state_tolerance = 1e-8
max_timesteps          = 120
model_time             = 0.0
max_timestep           = 1.0

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number
k                      = Constant((0,1))
X = x,y                = SpatialCoordinate(mesh)
h                      = sqrt(CellVolume(mesh))

# Temperature related constants:
delta_t                = Constant(5.0e-6) # Initial time-step
kappa                  = Constant(1.0)  # Thermal diffusivity
target_cfl_no          = 2.5
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
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_type': 'mumps',
    'mat_type': 'aij'
}

#########################################################################################################
################################## Geometry and Spatial Discretization: #################################
#########################################################################################################

# Set up function spaces - currently using the P2P1 element pair :
V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W    = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Q    = FunctionSpace(mesh, "CG", 1) # Temperature function space (scalar)

# Set up mixed function space and associated test functions:
Z    = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)
Y    = TestFunction(Q)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z    = Function(Z) # a field over the mixed function space Z.
u, p = split(z) # can we nicely name mixed function space fields?

# time step function
ts_func = Function(Q)

def compute_timestep(u):
    # A function to compute the timestep, based upon the CFL criterion
    ts_func.interpolate( h / sqrt(dot(u,u)))
    ts_min = ts_func.dat.data.min()
    ts_min = mesh.comm.allreduce(ts_min, MPI.MIN)
    return min(ts_min*target_cfl_no,max_timestep)

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################
final_state = Function(Q, name="ReferenceTemperature")
final_state_file = DumbCheckpoint("./FinalState", mode=FILE_READ)
final_state_file.load(final_state, name='Temperature')
final_state_file.close()

T_ic = Function(Q, name="InitialCondition")
T_checkpoint = DumbCheckpoint(basename='FinalState', mode=FILE_READ)
T_checkpoint.load(T_ic, name='Temperature')
T_checkpoint.close()

control = Control(T_ic)

# Set up temperature field and initialise based upon coordinates:
T_old   = Function(Q, name="OldTemperature")
T_old.assign(T_ic)

T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

#########################################################################################################
############################################ Setup Equations ############################################
#########################################################################################################

### Initially deal with Stokes equations ###

# Equation in weak (ufl) form - note that continuity equation is added here - need to better understand why:
# Set up in residual form to ensure non-linear solvers are used.
mu_field    = Function(W, name="Viscosity")
mu          = 1.0*exp(-ln(10)*T_new) # Variable viscosity
tau         = mu * (grad(u)+transpose(grad(u))) # Strain-rate tensor:
tauii       = sqrt((tau[0,0]**2 + tau[0,1]**2 + tau[0,1]**2 + tau[1,1]**2) / 2.)
F_stokes    = inner(grad(N), tau) * dx + dot(N,grad(p)) * dx - (dot(N,k)*Ra*T_theta) * dx
F_stokes   += dot(grad(M),u) * dx # Continuity equation

# Set up boundary conditions for Stokes: We first need to extract the velocity field from the mixed
# function space Z. This is done using Z.sub(0). We subsequently need to extract the x and y components
# of velocity. This is done using an additional .sub(0) and .sub(1), respectively. Note that the final arguments
# here are the physical boundary ID's.
u_top = Function(V, name="Velocity_top") # This is the velocity for the top boundary, that will be updated each time-step
bcv_fs_t      = DirichletBC(Z.sub(0), u_top, (top_id))
bcv_fs_b      = DirichletBC(Z.sub(0).sub(1), 0.0, (bottom_id))

# Generating nullspaces 
c0V = Function(V)
c0V.interpolate(Constant([1., 0.]))

V_nullspace = VectorSpaceBasis([c0V])
V_nullspace.orthonormalize()
Z_nullspace = MixedVectorSpaceBasis(Z, [V_nullspace, Z.sub(1)])

### Next deal with Temperature advection-diffusion equation: ###
delta_x      = sqrt(CellVolume(mesh))
Y_SUPG       = Y + dot(u,grad(Y)) * (delta_x / (2*sqrt(dot(u,u))))
F_energy     = Y_SUPG * ((T_new - T_old) / delta_t) * dx + Y_SUPG*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx
bct_base     = DirichletBC(Q, 1.0, bottom_id)
bct_top      = DirichletBC(Q, 0.0, top_id)

# Write output files in VTK format:
u, p = z.split() # Do this first to extract individual velocity and pressure fields:
u._name = 'Velocity' # This doesn't seem to work!
p._name = 'Pressure'

# Setting up the stokes solver
# Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
# and the RHS == 0.
stokes_prob = vs.NonlinearVariationalProblem(F_stokes, z, bcs=[bcv_fs_b, bcv_fs_t])
stokes_solver = vs.NonlinearVariationalSolver(stokes_prob,
                                              solver_parameters=stokes_solver_parameters,
                                              appctx={"mu": mu})

# Setting up the energy solver 
energy_prob = vs.NonlinearVariationalProblem(F_energy, T_new,
                                             bcs=[bct_base,bct_top])
energy_solver = vs.NonlinearVariationalSolver(energy_prob, solver_parameters=temperature_solver_parameters)



# Starting O file for storing Temperatures
checkpoint_boundary = DumbCheckpoint("./VelocityCheckpoint",mode=FILE_READ)

# Now perform the time loop:
for timestep in range(0,max_timesteps):
    log(f"Timestep Number: {timestep}, Time: {model_time}, delta_t: {float(delta_t)}")

    # Loading top boundary condition 
    checkpoint_boundary.set_timestep(t=model_time, idx=timestep)
    checkpoint_boundary.load(u_top, name='function_5[0]')

    # solve stokes
    stokes_solver.solve()

    # Let's keep delta_t constant 
    #delta_t.assign(compute_timestep(u))

    # solve temperature 
    energy_solver.solve()

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Update model_time
    model_time += float(delta_t)

checkpoint_boundary.close()



# define the objective functional
functional = assemble(0.5 * (T_new - final_state)**2 * dx)


# Below are callbacks allowing us to access various field information (accessed through reducedfunctional).
class OptimisationOutputCallbackPost:
    def __init__(self):
        self.iter_idx = 0
        self.opt_file             = File('visual/opt_file.pvd') 
        self.grad_copy            = Function(Q, name="Gradient")
        self.T_ic_copy            = Function(Q, name="InitTemperature")
        self.T_tc_copy            = Function(Q, name="FinTemperature")
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
        #
        # output current final state velocity and pressure
        self.z_copy.assign(z.block_variable.checkpoint)
        self.u_copy.assign(self.z_copy.split()[0])
        self.p_copy.assign(self.z_copy.split()[1])
        
        #  Write out the fields
        self.opt_file.write(self.T_ic_copy, self.T_tc_copy, self.u_copy, self.p_copy, self.grad_copy)
        func_val = assemble((self.T_tc_copy-final_state)**2 * dx) 
        #reg_val  = assemble(inner(grad(self.T_ic_copy-T_mean), grad(self.T_ic_copy - T_mean)) * dx) 
        grad_val = assemble((self.grad_copy)**2 * dx)

        log(f"# Der {self.iter_idx}, ||Misfit|| {func_val}, ||Grad|| {grad_val}")
        log('Derivative calculation', self.iter_idx)
        self.iter_idx += 1

class ForwardCallbackPost:
   def __init__(self):
      self.fwd_idx = 0 
   def __call__(self, func_value, controls):
      self.fwd_idx +=1 
      log('# fwd {}, ||Misfit|| {}'.format(self.fwd_idx, func_value))

# Initiate classes for the callbacks
local_cb_post = OptimisationOutputCallbackPost()
eval_cb_post = ForwardCallbackPost()

# Defining the object for pyadjoint
reduced_functional = ReducedFunctional(functional, control, eval_cb_post=eval_cb_post, derivative_cb_post=local_cb_post)

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



# This is the classic way
params = {
        'General': {
                'Print Verbosity':1,
                'Secant': {'Type': 'Limited-Memory BFGS', 'Maximum Storage': 5}, 
                    },
        'Step': {
           'Type': 'Line Search',
           'Line Search': {
                'Descent Method': {'Type': 'Quasi-Newton Method'},
                'Line-Search Method': {'Type': 'Cubic Interpolation'}, 
                #'Function Evaluation Limit': 5,
                #'Sufficient Decrease Tolerance': 1e-2,
                #'Use Previous Step Length as Initial Guess': True,
                            }
                },
        'Status Test': {
            'Gradient Tolerance': 1e-12,
            'Iteration Limit': 100,
                        }
        }

rol_solver = ROLSolver(minp, params)
params = ROL.ParameterList(params, "Parameters")
rol_secant = InitHessian(5) # maximum storage
#rol_step = ROL.TrustRegionStep(rol_secant, params)
rol_step = ROL.LineSearchStep(params, rol_secant)
rol_status = ROL.StatusTest(params)
rol_algorithm = ROL.Algorithm(rol_step, rol_status)

with stop_annotating():
    rol_algorithm.run(
        rol_solver.rolvector,
        rol_solver.rolobjective,
        rol_solver.bounds # only if we have a bounded problem
            )
    optimal_ic = rol_solver.problem.reduced_functional.controls.delist(rol_solver.rolvector.dat)
    # Save the optimal temperature_ic field 
    ckpt_T_ic = DumbCheckpoint("T_ic_optimal",\
            single_file=True, mode=FILE_CREATE,\
                               comm=mesh.comm)
    ckpt_T_ic.store(optimal_ic)
    ckpt_T_ic.close()





