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
import ROL as ROL
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
steady_state_tolerance = 1e-7
max_num_timesteps      = 1
target_cfl_no          = 2.5
max_timestep           = 1.00

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
V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W    = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Q    = FunctionSpace(mesh, "CG", 1) # Temperature function space (scalar)

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

# Final state, which will be used as reference for minimization, loaded from a file 
final_state = Function(Q, name='RefTemperature')
final_state.assign(0.1)

# Initial condition
T_ic   = Function(Q, name="T_IC")
# Let's start with the final condition
T_ic.assign(0.6)

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

# Stokes in weak form 
F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
F_stokes += - (dot(N,yhat)*Ra*T_theta) * dx 
F_stokes += - div(u)* M * dx

# Setting No-slip BC for top and bottom
bcu_topbase     = DirichletBC(Z.sub(0), 0.0, (top_id, bottom_id))
bcu_rightleft   = DirichletBC(Z.sub(0), 0.0, (left_id, right_id))

# Pressure nullspace                                                                                                                                                                                                                                                                      
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

### Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# For some reason this only works here!
u, p    = z.split() # Do this first to extract individual velocity, pressure and lagrange multplier fields:
u.rename('Velocity') 
p.rename('Pressure')

# A simulation time to track how far we are
simu_time = 0.0

z_tri = TrialFunction(Z)
F_stokes_lin = replace(F_stokes, {z: z_tri})
a, L = lhs(F_stokes_lin), rhs(F_stokes_lin)
stokes_problem = LinearVariationalProblem(a, L, z, constant_jacobian=True, bcs=[bcu_topbase, bcu_rightleft])
stokes_solver  = LinearVariationalSolver(stokes_problem, solver_parameters=solver_parameters, nullspace=p_nullspace, transpose_nullspace=p_nullspace)

q_tri = TrialFunction(Q)
F_energy_lin = replace(F_energy, {T_new:q_tri})
a_energy, L_energy = lhs(F_energy_lin), rhs(F_energy_lin)
energy_problem = LinearVariationalProblem(a_energy, L_energy, T_new, constant_jacobian=False)
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

    # Updating Temperature
    log("Timestep Number: ", timestep, " Timestep: ", float(delta_t))

## Initialise functional
functional = assemble(0.5*(T_new - final_state)**2 * dx)


# Defining the object for pyadjoint
reduced_functional = ReducedFunctional(functional, control)

### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional)

# This is the classic way
params = {
        'General': {
              'Print Verbosity': 1,
              'Secant': {'Type': 'Limited-Memory BFGS', 'Maximum Storage': 10},
                    },
        'Step': {
           'Type': 'Line Search',
           'Line Search': {
                'Descent Method': {'Type': 'Quasi-Newton Method'},
                'Line-Search Method': {
                                'Type':  "Brent's", #'Cubic Interpolation ''Backtracking',#'Bisection',
                                'Backtracking Rate': 0.5,
                                'Bracketing Tolerance': 1.e-1,
                                'Bisection': {
                                    'Tolerance': 1e-1,
                                    'Iteration Limit': 20,
                                            },
                                        },
                                "Brent's": {
                                    'Tolerance': 1e0,
                                    'Iteration Limit': 3,
                                    'Run Test Upon Initialization': False,
                                            },
                'Curvature Condition': {
                                'Type': 'Strong Wolfe Conditions',
                                'General Parameter': 0.9,
                                'Generalized Wolfe Parameter': 0.6,
                                        },
                'Function Evaluation Limit': 20,
                'Sufficient Decrease Tolerance': 1e-1,
                'Use Previous Step Length as Initial Guess': False,
                            },
                },
        'Status Test': {
            'Gradient Tolerance': 0,
            'Iteration Limit': 1,
                        }
        }


with stop_annotating():
    # set up ROL problem
    rol_solver = ROLSolver(minp, params, inner_product="L2")
    sol = rol_solver.solve()


