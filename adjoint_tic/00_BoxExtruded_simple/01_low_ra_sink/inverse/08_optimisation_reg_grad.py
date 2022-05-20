"""
    Adjoint Reconstruction - Using the classic way inputing parameters, instead of definiting methods for ROL.Algirithm() 
"""
from re import I
import sys
sys.path.append('../../')

from firedrake import *
from mpi4py import MPI
from firedrake.petsc import PETSc
from lib_averaging import LayerAveraging
from firedrake_adjoint import *
from pyadjoint import MinimizationProblem, ROLSolver
from pyadjoint.tape import no_annotations, Tape, set_working_tape
import ROL 
import time

alpha = 0.01
simu_name = f"08_reg_grad_{alpha}"

# Quadrature degree:
dx = dx(degree=6)

# logging.set_log_level(1)
# logging.set_level(1)

with CheckpointFile("../../Final_State.h5", "r") as T_ref_checkpoint:
    mesh = T_ref_checkpoint.load_mesh("firedrake_default_extruded")
    final_state = T_ref_checkpoint.load_function(mesh, "Temperature")
    T_average = T_ref_checkpoint.load_function(mesh, "OneDimTemperature")
    final_state.rename("final_state")

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'
left_id, right_id = 1, 2

# spatial coordinates
X = x, y = SpatialCoordinate(mesh)

# setting up vertical direction
y_abs = sqrt(y**2)
yhat = as_vector((0, y)) / y_abs

# Global Constants:
max_num_timesteps = 70

# Stokes related constants:
Ra = Constant(1e6)  # Rayleigh Number

# Below are callbacks relating to the adjoint solutions
# (accessed through solve).
# Not sure what the best place would be to initiate working tape!
tape = get_working_tape()

# Temperature related constants:
delta_t = Constant(5e-6)  # Time-step
kappa = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson
# scheme where theta_ts = 0.5.
theta_ts = 0.5


# Print function to ensure log output is only written
# on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args)


# Geometry and Spatial Discretization:
# Set up function spaces - currently using the P2P1 element pair :
V = VectorFunctionSpace(mesh, "CG", 2)  # Velocity function space (vector)
W = FunctionSpace(mesh, "CG", 1)  # Pressure function space (scalar)
Q = FunctionSpace(mesh, "CG", 2)  # Temperature function space (scalar)

# Set up mixed function space and associated test functions:
Z = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)
Y = TestFunction(Q)

# Set up fields on these function spaces - split into each component
# so that they are easily accessible:
z = Function(Z)  # a field over the mixed function space Z.
u, p = split(z)     # can we nicely name mixed function space fields?

# Setting free-slip BC for bottom, and sides
u_surface = Function(V, name="BoundaryCondition")

# T advection diffusion equation Prerequisites
# Initial condition, let's start with the final condition
T_ic = Function(Q, name="T_IC")
T_ic.interpolate(final_state)

# Set up temperature field and initialise it with present day
T_old = Function(Q, name="OldTemperature")
T_old.assign(T_ic)

T_new = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

# Setup Equations
# viscosiy
mu = Constant(1.0)  # Constant viscosity


# deviatoric stresses
def tau(u):
    return mu * (grad(u) + transpose(grad(u)))


# Stokes in weak form
F_stokes = inner(grad(N), tau(u)) * dx - div(N)*p * dx
F_stokes += - (dot(N, yhat)*Ra*T_theta) * dx
F_stokes += - div(u) * M * dx

# Setting free-slip BC for top and bottom
bcu_topbase = DirichletBC(Z.sub(0).sub(1), 0.0, (bottom_id, top_id))
bcu_rightleft = DirichletBC(Z.sub(0).sub(0), 0.0, (left_id, right_id))
all_u_bounds = [bcu_topbase, bcu_rightleft]

# Setting Dirihlet BC for top and bottom of the temperature field
bct_top = DirichletBC(Q, 0.0, (top_id))
bct_base = DirichletBC(Q, 1.0, (bottom_id))

# Pressure nullspace
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0),
                                     VectorSpaceBasis(constant=True)]
                                    )

# Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx\
    + Y*dot(u,grad(T_theta)) * dx\
    + dot(grad(Y),kappa*grad(T_theta)) * dx

# Setup problem and solver objects so we can reuse (cache) solver setup
stokes_problem = NonlinearVariationalProblem(F_stokes, z, bcs=all_u_bounds)
stokes_solver = NonlinearVariationalSolver(stokes_problem,
                                           nullspace=p_nullspace,
                                           transpose_nullspace=p_nullspace,
                                           )

energy_problem = NonlinearVariationalProblem(F_energy, T_new,
                                             bcs=[bct_base, bct_top]
                                             )
energy_solver = NonlinearVariationalSolver(energy_problem)

# Setting adjoint and forward callbacks, and control parameter
control = Control(T_ic)

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems,
    # where everything is on the LHS (as above)
    # and the RHS == 0.
    stokes_solver.solve()

    # Temperature system:
    energy_solver.solve()

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

# Initialise functional
functional = assemble(0.5*(T_new - final_state)**2 * dx)/assemble(0.5*(final_state)**2 * dx)
regularisation = assemble(0.5*(dot(grad(T_ic - T_average), grad(T_ic - T_average))) * dx)\
                /assemble(0.5*(dot(grad(T_average), grad(T_average))) * dx)

class myReducedFunctional(ReducedFunctional):
    def __init__(self, functional, controls, **kwargs):
        super().__init__(functional, control, **kwargs)
        self.fwd_cntr = 0
        self.adj_cntr = 0

    def __call__(self, values):
        pre_time = time.time()
        val = super().__call__(values)
        self.fwd_cntr += 1
        log(f"\tFWD # {self.fwd_cntr} call took {time.time() - pre_time}")
        return val

    def derivative(self):
        pre_time = time.time()
        deriv = super().derivative(options={})
        self.adj_cntr += 1
        log(f"\tADJ # {self.adj_cntr} call took {time.time() - pre_time}")
        return deriv


# Defining the object for pyadjoint
reduced_functional = myReducedFunctional(functional + alpha * regularisation,
                                         control
                                        )

# Set up bounds, which will later be used to
# enforce boundary conditions in inversion:
T_lb = Function(Q, name="LB_Temperature")
T_ub = Function(Q, name="UB_Temperature")
T_lb.assign(0.0)
T_ub.assign(1.0)

# Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))


class myStatusTest(ROL.StatusTest):
    def __init__(self, params, vector):
        super().__init__(params)

        # Keep track of the vector that is being passed to StatusCheck
        self.vector = vector
        self.T_copy = Function(Q, name="Temperature")
        self.opt_t_final_file = File(f'visual_{simu_name}/opt_temperature_fin.pvd')
        self.opt_t_init_file = File(f'visual_{simu_name}/opt_temperature_int.pvd')

        self.solution_checkpoint = CheckpointFile(f"solution_{simu_name}.h5", mode="w")
        self.solution_checkpoint.save_mesh(mesh)

        # loading the true answer for comparisons 
        with CheckpointFile("../../Initial_State.h5", 'r') as true_initial_state_file:
            _ = true_initial_state_file.load_mesh("firedrake_default_extruded")
            self.T_ic_true = true_initial_state_file.load_function(mesh, "Temperature")
            self.T_ic_true.rename("ref_initial_state")
        self.my_idx = 0

    @no_annotations
    def check(self, status):

        log("\tFinal misfit: {}".format(assemble(0.5*(T_new.block_variable.checkpoint - final_state)**2 * dx)))
        log("\tInitial misfit : {}".format(assemble(0.5*(self.vector.dat[0] - self.T_ic_true)**2 * dx)))
        # log("regularisation: {}".format(assemble(0.5*beta*dot(grad(self.vector.dat[0]), grad(self.vector.dat[0])) * dx)))

        # Writing out final condition
        self.T_copy.assign(T_new.block_variable.checkpoint)
        self.opt_t_final_file.write(self.T_copy, final_state)

        # Writing out initial condition
        self.T_copy.assign(self.vector.dat[0])
        self.opt_t_init_file.write(self.T_copy, self.T_ic_true)

        # Write out the solution
        self.solution_checkpoint.save_function(self.T_copy, idx=self.my_idx)
        self.my_idx += 1

        return ROL.StatusTest.check(self, status)


# This is the classic way
params = {
        'General': {
              'Print Verbosity': 1,
              'Output Level': 3,
              'Krylov': {   # These are needed for our
                            # solution of the hessian I guess
                    "Iteration Limit": 10,
                    "Absolute Tolerance": 1e-4,
                    "Relative Tolerance": 1e-2,
                    },
              'Secant': {'Type': 'Limited-Memory BFGS',
                         'Maximum Storage': 20,
                         'Use as Hessian': True,
                         "Barzilai-Borwein": 1},
                    },
        'Step': {
           'Type': 'Trust Region',  # 'Line Search',
           'Trust Region': {
                "Lin-More":     {
                    "Maximum Number of Minor Iterations": 10,
                    "Sufficient Decrease Parameter":      1e-2,
                    "Relative Tolerance Exponent":        1.0,
                    "Cauchy Point": {
                        "Maximum Number of Reduction Steps": 10,
                        "Maximum Number of Expansion Steps": 10,
                        "Initial Step Size":                 1.0,
                        "Normalize Initial Step Size":       False,
                        "Reduction Rate":                    0.1,
                        "Expansion Rate":                    10.0,
                        "Decrease Tolerance":                1e-8,
                                    },
                        "Projected Search": {
                                "Backtracking Rate": 0.5,
                                "Maximum Number of Steps": 20,
                                            },
                                },
                #  Subproblem Model could be "Kelley-Sachs",
                "Subproblem Model":                     "Lin-More",
                "Initial Radius":                       0.005,
                "Maximum Radius":                       1e20,
                "Step Acceptance Threshold":            0.05,
                "Radius Shrinking Threshold":           0.05,
                "Radius Growing Threshold":             0.9,
                "Radius Shrinking Rate (Negative rho)": 0.0625,
                "Radius Shrinking Rate (Positive rho)": 0.25,
                "Radius Growing Rate":                  2.5,
                "Sufficient Decrease Parameter":        1.e-2,
                "Safeguard Size":                       100,
                            },
                },
        'Status Test': {
            'Gradient Tolerance': 0,
            'Iteration Limit': 50,
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
