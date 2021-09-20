"""
    Adjoint Reconstruction - Using the classic way inputing parameters, instead of definiting methods for ROL.Algirithm() 
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
max_num_timesteps      = 40
target_cfl_no          = 2.5
max_timestep           = 1.00

# Stokes related constants:
Ra                     = Constant(1e8)   # Rayleigh Number

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
final_state_file = DumbCheckpoint("../../final_state", mode=FILE_READ)
final_state_file.load(final_state, 'Temperature')
final_state_file.close()

T_mean = Constant(0.5)

# Initial condition
T_ic   = Function(Q, name="T_IC")
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
#bcu_topbase     = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id)) # bottom boundary is always free-slip
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
boundary_velocity_fi = DumbCheckpoint("../../velocity", mode=FILE_READ)

# Documenting the first forward run
fwd_run_fi = File('visual/first_fwd_run.pvd')

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
    solve(F_energy==0, T_new, solver_parameters=stokes_solver_parameters)

    fwd_run_fi.write(T_new, u, p)

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log("Timestep Number: ", timestep, " Timestep: ", str('%10.4e' %simu_time))

boundary_velocity_fi.close()



## Initialise functional
#beta_m = 1.0/assemble((final_state - T_mean)**2 * dx)
#misfit = assemble(float(beta_m) * (T_new - final_state)**2 * dx)
misfit = assemble((T_new - final_state)**2 * dx)
#
#beta_r = 1e-3/ assemble(inner(grad(final_state - T_mean), grad(final_state - T_mean)) * dx)
#reg    = assemble(float(beta_r) * inner(grad(T_ic - T_mean), grad(T_ic - T_mean)) * dx)

# functional is a combination of the misfit and the regularising term
#functional = misfit + reg
functional = misfit 

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
        reg_val  = assemble(inner(grad(self.T_ic_copy-T_mean), grad(self.T_ic_copy - T_mean)) * dx) 
        grad_val = assemble((self.grad_copy)**2 * dx)

        log('# Der {}, ||Misfit|| {}, Non-existing ||Regu|| {}, ||Grad|| {}'.format(self.iter_idx, func_val, reg_val, grad_val))
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
            'Iteration Limit': 500,
                        }
        }


with stop_annotating():    
    # set up ROL problem
    rol_solver = ROLSolver(minp, params, inner_product="L2")
    sol = rol_solver.solve()

    # Save the optimal temperature_ic field 
    ckpt_T_ic = DumbCheckpoint("T_ic_optimal",\
            single_file=True, mode=FILE_CREATE,\
                               comm=mesh.comm)
    ckpt_T_ic.store(sol)
    ckpt_T_ic.close()

