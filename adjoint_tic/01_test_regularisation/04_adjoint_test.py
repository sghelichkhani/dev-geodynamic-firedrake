"""
    Testing the derivative of regularisation
"""
import sys
from firedrake import *
from mpi4py import MPI
from firedrake.petsc import PETSc
from lib_averaging import LayerAveraging
from firedrake_adjoint import *
from pyadjoint import MinimizationProblem, ROLSolver, ROLVector
from pyadjoint.tape import no_annotations, Tape, set_working_tape
import ROL 
import time
import numpy as np

# Quadrature degree:
dx = dx(degree=6)


with CheckpointFile("./Final_State.h5", "r") as T_ref_checkpoint:
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
max_num_timesteps = 0  
simu_time = 0.0

# Stokes related constants:
Ra = Constant(1e6)  # Rayleigh Number

# Below are callbacks relating to the adjoint solutions
# (accessed through solve).
# Not sure what the best place would be to initiate working tape!
tape = get_working_tape()

# Temperature related constants:
delta_t = Constant(4.0e-6)  # Time-step
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

# Setting Dirihlet BC for top and bottom of the temperature field
bct_top = DirichletBC(Q, 0.0, (top_id))
bct_base = DirichletBC(Q, 1.0, (bottom_id))
all_t_bounds = [bct_top, bct_base]

# T advection diffusion equation Prerequisites
# Initial condition, let's start with the final condition
T_ic = Function(Q, name="T_IC")
T_ic.interpolate(final_state)

# Set up temperature field and initialise it with present day
T_old = Function(Q, name="OldTemperature")
T_old.assign(T_ic)

T_average_Q2 = Function(Q, name="Q2field")
T_average_Q2.project(T_average)

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


# Pressure nullspace
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0),
                                     VectorSpaceBasis(constant=True)]
                                    )
# Generating near_nullspaces for GAMG:
z_rotV = Function(V).interpolate(as_vector((-X[1], X[0])))
nns_x = Function(V).interpolate(Constant([1., 0.]))
nns_y = Function(V).interpolate(Constant([0., 1.]))
V_near_nullspace = VectorSpaceBasis([nns_x, nns_y, z_rotV])
V_near_nullspace.orthonormalize()
Z_near_nullspace = MixedVectorSpaceBasis(Z, [V_near_nullspace, Z.sub(1)])

# Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx\
    + Y*dot(u,grad(T_theta)) * dx\
    + dot(grad(Y),kappa*grad(T_theta)) * dx

# Setup problem and solver objects so we can reuse (cache) solver setup
stokes_problem = NonlinearVariationalProblem(F_stokes, z, bcs=all_u_bounds)
stokes_solver = NonlinearVariationalSolver(stokes_problem, 
     nullspace=p_nullspace, transpose_nullspace=p_nullspace,
)

energy_problem = NonlinearVariationalProblem(F_energy, T_new, bcs=[bct_base, bct_top])
energy_solver = NonlinearVariationalSolver(energy_problem)


# Setting adjoint and forward callbacks, and control parameter
control = Control(T_ic)

# Make sure boundary conditions are applied
solve(Y*TrialFunction(Q)*dx == Y*T_ic*dx, T_ic, bcs=[bct_base, bct_top])

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

# Compute the analytical expression for the gradient term
analyticalder = Function(Q, name="AnalyticalDer").interpolate(-div(grad(T_ic-T_average)))


# Initialise functional
functional = assemble(0.5*(T_new - final_state)**2 * dx)
regularisation = assemble(0.5*(inner(grad(T_ic-T_average), grad(T_ic-T_average))) * dx)


# Defining the object for pyadjoint, we use regulatisation seperately
reduced_functional = ReducedFunctional(regularisation, control)

reduced_functional([T_ic])
pyadjointder = reduced_functional.derivative(options={'riesz_representation':'L2'})
pyadjointder.rename("PyAdjointDer")

# Compute analytical term
myfi = File("RegDerivatives.pvd")
myfi.write(analyticalder, pyadjointder)

#
#Delta_temp  = Function(Q, name="Delta_Temperature")
#Delta_temp.dat.data[:] = np.random.random(Delta_temp.dat.data.shape)
#minconv = taylor_test(reduced_functional, T_ic, Delta_temp)
#log(minconv)


# Some tests to be done with these
#invmass = Function(Q, name="InverseMass")
#mass_form = Y*TrialFunction(Q)*dx
#
#mass_action_form = assemble(action(mass_form, Constant(1)))
#
#ls = LinearSolver(assemble(mass_form))
#ls.solve(invmass, interpolate(Constant(1.0), Q))
#


