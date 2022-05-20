"""
    Adjoint Reconstruction
    We try a tikhonov regularisation here
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
x_max, y_max = 1.0, 1.0
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

# setting up vertical direction
y_abs     = sqrt(y**2)
yhat  = as_vector((0,y)) / y_abs

# Global Constants:
max_num_timesteps      = 60

# Stokes related constants:
Ra                     = Constant(1e7)   # Rayleigh Number

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


# Final states the like of tomography, note that we also load the reference profile 
final_state = Function(Q, name='RefTemperature')

final_state_file = DumbCheckpoint("../../forward/State_100", mode=FILE_READ)
final_state_file.load(final_state, 'Temperature')
final_state_file.close()

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

# Setup problem and solver objects so we can reuse (cache) solver setup
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

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

## Initialise functional
alpha = 0.001
functional = assemble(0.5*(T_new - final_state)**2 * dx)
regularisation = assemble(alpha * 0.5*(T_ic - 0.5)**2 * dx)


class myReducedFunctional(ReducedFunctional):
    def __init__(self, functional, controls, **kwargs):
        super().__init__(functional, control, **kwargs)
        self.fwd_cntr = 0
        self.adj_cntr = 0

    def __call__(self, values):
        pre_time = time.time()
        val = super().__call__(values)
        self.fwd_cntr +=1
        log(f"\tFWD # {self.fwd_cntr} call took {time.time() - pre_time}")
        return val

    def derivative(self):
        pre_time = time.time()
        deriv = super().derivative(options={})
        self.adj_cntr +=1
        log(f"\tADJ # {self.adj_cntr} call took {time.time() - pre_time}")
        return deriv

# Defining the object for pyadjoint
reduced_functional = myReducedFunctional(functional+regularisation, control)

# Set up bounds, which will later be used to enforce boundary conditions in inversion:
T_lb     = Function(Q, name="LB_Temperature")
T_ub     = Function(Q, name="UB_Temperature")
T_lb.assign(0.0)
T_ub.assign(1.0)

### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))

class myStatusTest(ROL.StatusTest):
    def __init__(self, params, vector):
        super().__init__(params)
        #ROL.StatusTest.__init__(self, params)

        # Keep track of the vector that is being passed to StatusCheck
        self.vector = vector
        self.T_copy              = Function(Q, name="Temperature")
        self.opt_t_final_file    = File(str('visual_beta_%5.4f_%3.3i/opt_temperature_fin.pvd' %(alpha, max_num_timesteps)))
        self.opt_t_init_file     = File(str('visual_beta_%5.4f_%3.3i/opt_temperature_int.pvd' %(alpha, max_num_timesteps)))

        # A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
        self.T_ic_true              = Function(Q, name="InitTemperature_Ref")
        true_initial_state_file = DumbCheckpoint(str("../../forward/State_%3.3i" %(100-max_num_timesteps)) , mode=FILE_READ)
        true_initial_state_file.load(self.T_ic_true,    'Temperature')
        true_initial_state_file.close()

    @no_annotations
    def check(self, status):

        log("\tFinal misfit: {}".format(assemble(0.5*(T_new.block_variable.checkpoint - final_state)**2 * dx)/assemble(0.5*(final_state **2) * dx)))
        log("\tInitial misfit : {}".format(assemble(0.5*(self.vector.dat[0] - self.T_ic_true)**2 * dx)))
        log("\tRegularisation: {}".format(assemble(alpha * 0.5*(self.vector.dat[0] - 0.5)**2 * dx)))
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
              'Output Level': 3,
              'Krylov': {# These are needed for our solution of the hessian I guess
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
           'Type': 'Trust Region',#'Line Search',
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
                "Subproblem Model":                     "Lin-More",#"Kelley-Sachs",
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
    # Save the optimal temperature_ic field 
    ckpt_T_ic = DumbCheckpoint(str("T_ic_State_%3.3i" %(max_num_timesteps)),\
            single_file=True, mode=FILE_CREATE,\
                               comm=mesh.comm)
    ckpt_T_ic.store(optimal_ic)
    ckpt_T_ic.close()


