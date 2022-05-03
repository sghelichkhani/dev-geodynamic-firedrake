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
import time

#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

#logging.set_log_level(1)
#logging.set_level(1)

# Geometric Constants:
y_max = 1.0
x_max = 1.0

#  how many intervals along x/y directions 
disc_n = 100

# and Interval mesh of unit size 
mesh1d = IntervalMesh(disc_n, length_or_left=0.0, right=x_max) 
# extruding the base mesh "mesh1d" in the third dimension
mesh = ExtrudedMesh(mesh=mesh1d, layers=disc_n, layer_height=y_max/disc_n, extrusion_type='uniform', kernel=None, gdim=None)

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'
left_id, right_id = 1, 2

# spatial coordinates
X  = x, y = SpatialCoordinate(mesh)

# a measure of cell size
h	  = sqrt(CellVolume(mesh))

# setting up vertical direction
y_abs     = sqrt(y**2)
yhat  = as_vector((0,y)) / y_abs

# Global Constants:
max_num_timesteps      = 20

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number

# Below are callbacks relating to the adjoint solutions (accessed through solve).
# Not sure what the best place would be to initiate working tape!
tape = get_working_tape()

# Temperature related constants:
delta_t                = Constant(5e-6) # Time-step
kappa                  = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson scheme where theta_ts = 0.5:
theta_ts               = 0.5

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args) 

# Stokes Equation Solver Parameters:
solver_parameters = {
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
V       = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W       = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Q       = FunctionSpace(mesh, "CG", 2) # Temperature function space (scalar)
Qlayer  = FunctionSpace(mesh1d, "CG", 2)

# Set up mixed function space and associated test functions:
Z       = MixedFunctionSpace([V, W])
N, M    = TestFunctions(Z)
Y       = TestFunction(Q)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z    = Function(Z)  # a field over the mixed function space Z.
u, p = split(z)     # can we nicely name mixed function space fields?

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Reference Initial State
# This will be used for "True Misfit" calculations
true_initial_state = Function(Q, name='TrueInitialState')

# Final states the like of tomography, note that we also load the reference profile 
final_state = Function(Q, name='RefTemperature')
T_averef    = Function(Q, name='RefProfile')

final_state_file = DumbCheckpoint("../../forward/State_100", mode=FILE_READ)
final_state_file.load(final_state, 'Temperature')
final_state_file.load(T_averef,    'LayerAveragedTemp')
final_state_file.close()

# Radoially averaged temperature field
T_ave       = Function(Q, name='TemperatureProfile')

# Initial condition, let's start with the final condition
T_ic   = Function(Q, name="T_IC")
T_ic.project(final_state)

# Set up temperature field and initialise it with present day 
T_old    = Function(Q, name="OldTemperature")
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

# Stokes in weak form 
F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
F_stokes += - (dot(N,yhat)*Ra*T_theta) * dx 
F_stokes += - div(u)* M * dx

# Setting No-slip BC for top and bottom
bcu_topbase     = DirichletBC(Z.sub(0), 0.0, (top_id, bottom_id))
bcu_rightleft   = DirichletBC(Z.sub(0), 0.0, (left_id, right_id))

# Pressure nullspace
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

# Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# Setting Dirichlet temperature boundary condition
bcT_top     = DirichletBC(Q, 0.0, (top_id))
bcT_base    = DirichletBC(Q, 1.0, (bottom_id))

# For some reason this only works here!
u, p    = z.split() # Do this first to extract individual velocity, pressure and lagrange multplier fields:
u.rename('Velocity') 
p.rename('Pressure')

# A simulation time to track how far we are
simu_time = 0.0

# Setup problem and solver objects so we can reuse (cache) solver setup
z_tri = TrialFunction(Z)
F_stokes_lin = replace(F_stokes, {z: z_tri})
a, L = lhs(F_stokes_lin), rhs(F_stokes_lin)
stokes_problem = LinearVariationalProblem(a, L, z, constant_jacobian=True, bcs=[bcu_topbase, bcu_rightleft])
stokes_solver  = LinearVariationalSolver(stokes_problem, solver_parameters=solver_parameters, nullspace=p_nullspace, transpose_nullspace=p_nullspace)

q_tri = TrialFunction(Q)
F_energy_lin = replace(F_energy, {T_new:q_tri})
a_energy, L_energy = lhs(F_energy_lin), rhs(F_energy_lin)
energy_problem = LinearVariationalProblem(a_energy, L_energy, T_new, constant_jacobian=False, bcs=[bcT_top, bcT_base])
energy_solver  = LinearVariationalSolver(energy_problem, solver_parameters=solver_parameters)

# Setting adjoint and forward callbacks, and control parameter
control = Control(T_ic)

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    stokes_solver.solve()

    # Temperature system:
    energy_solver.solve()

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

## Initialise functional
functional = assemble(0.5*(T_new - final_state)**2 * dx)

# Below are callbacks allowing us to access various field information (accessed through reducedfunctional).
class OptimisationOutputCallbackPost:
    def __init__(self):
        self.iter_idx = 0
    def __call__(self, cb_functional, dj, controls):
        log(f'# Der {self.iter_idx}')
        self.iter_idx += 1

class ForwardCallbackPost:
   def __init__(self):
      self.fwd_idx = 0 
   def __call__(self, func_value, controls):
      self.fwd_idx +=1 
      log(f'# fwd {self.fwd_idx}, ||Misfit|| {func_value}')

class myReducedFunctional(ReducedFunctional):
    def __call__(self, values):
        pre_time = time.time()
        val = super().__call__(values)
        log(f"\tFWD call took {time.time() - pre_time}")
        return val 
    def derivative(self):
        pre_time = time.time()
        deriv = super().derivative(options={})
        log(f"\tADJ call took {time.time() - pre_time}")
        return deriv

# Initiate classes for the callbacks
local_cb_post = OptimisationOutputCallbackPost()
eval_cb_post = ForwardCallbackPost()

# Defining the object for pyadjoint
reduced_functional = myReducedFunctional(functional, control, eval_cb_post=eval_cb_post, derivative_cb_post=local_cb_post)

# Set up bounds, which will later be used to enforce boundary conditions in inversion:
T_lb     = Function(Q, name="LB_Temperature")
T_ub     = Function(Q, name="UB_Temperature")
T_lb.assign(0.0)
T_ub.assign(1.0)

### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))

class myStatusTest(ROL.StatusTest):
    def __init__(self, params, vector):
        #super().__init__(self, params)
        ROL.StatusTest.__init__(self, params)
        # times will be used for measuring performance
        self.time_zero = time.time()
        self.time_past = time.time()

        # Keep track of the vector that is being passed to StatusCheck
        self.vector = vector
        self.T_copy              = Function(Q, name="Temperature")
        self.opt_t_final_file    = File('visual_short/opt_temperature_fin.pvd')
        self.opt_t_init_file     = File('visual_short/opt_temperature_int.pvd')

        # A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
        self.T_ic_true              = Function(Q, name="InitTemperature_Ref")
        true_initial_state_file = DumbCheckpoint("../../forward/State_40", mode=FILE_READ)
        true_initial_state_file.load(self.T_ic_true,    'Temperature')
        true_initial_state_file.close()

    @no_annotations
    def check(self, status):
        import resource
        log(f"\tru_maxrss {(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss + resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1e6}.")
        log(f"\tWall time {time.time()-self.time_zero}.")
        log(f"\tWall time for iteration {time.time()-self.time_past}.")
        self.time_past = time.time()

        log("\tFinal misfit: {}".format(assemble(0.5*(T_new.block_variable.checkpoint - final_state)**2 * dx)))
        log("\tInitial misfit : {}".format(assemble(0.5*(self.vector.dat[0] - self.T_ic_true)**2 * dx)))
        #log("regularisation: {}".format(assemble(0.5*beta*dot(grad(self.vector.dat[0]), grad(self.vector.dat[0])) * dx)))

        # Writing out final condition
        self.T_copy.assign(T_new.block_variable.checkpoint)
        self.opt_t_final_file.write(self.T_copy, final_state)

        ## Writing out initial condition
        self.T_copy.assign(self.vector.dat[0])
        self.opt_t_init_file.write(self.T_copy, self.T_ic_true)

        return ROL.StatusTest.check(self, status)


# This is the classic way                    
params = {                                   
        'General': {                   
              'Print Verbosity': 1,         
              'Output Level': 1,
              'Secant': {'Type': 'Limited-Memory BFGS',
                         'Maximum Storage': 20, 
                         'Use as Hessian': True,
                         "Barzilai-Borwein": 1,
                        },
                    },
        'Step': {
           'Type': 'Trust Region',  #'Line Search',
           'Trust Region': {
                "Subproblem Solver":                    "Truncated CG", #"Truncated CG",
                "Subproblem Model":                     "Lin-More", #"Kelley-Sachs",
                "Initial Radius":                       10.0,
                "Maximum Radius":                       1.e20,
                "Step Acceptance Threshold":            0.05,
                "Radius Shrinking Threshold":           0.05,
                "Radius Growing Threshold":             0.9,
                "Radius Shrinking Rate (Negative rho)": 0.0625,
                "Radius Shrinking Rate (Positive rho)": 0.25,
                "Radius Growing Rate":                  2.5,
                "Sufficient Decrease Parameter":        1.e-2,
                "Safeguard Size":                       1.e8,
                        },
                },  
        'Status Test': {
            'Gradient Tolerance': 0,
            'Iteration Limit': 5,
                        }
        }


rol_solver = ROLSolver(minp, params)
params = ROL.ParameterList(params, "Parameters")
status_test = myStatusTest(params, rol_solver.rolvector)

secant = ROL.InitBFGS(20)
rol_algorithm = ROL.LinMoreAlgorithm(params, secant)
rol_algorithm.setStatusTest(status_test, False)

with stop_annotating():
    rol_algorithm.run(
        rol_solver.rolvector,
        rol_solver.rolobjective,
        rol_solver.bounds # only if we have a bounded problem
            )
    optimal_ic = rol_solver.problem.reduced_functional.controls.delist(rol_solver.rolvector.dat)
    # Save the optimal temperature_ic field 
    ckpt_T_ic = DumbCheckpoint("T_ic_optimal_mid",\
            single_file=True, mode=FILE_CREATE,\
                               comm=mesh.comm)
    ckpt_T_ic.store(optimal_ic)
    ckpt_T_ic.close()


