""""
    Generation of forward problem - for an adjoint simulation
    The test problem is no heat-flow at boundaries, with free-slip on all surfaces
"""
import matplotlib.pyplot as plt
from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

# This will tell us how thick the 
thickness_val = 3.

# Some important constants etc
# Quadrature degree:
dx = dx(degree=6)

#
mesh = Mesh("t1_transfinite.msh")

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 3, 1 
left_id, right_id = 4, 2 

# spatial coordinates
X = x, y = SpatialCoordinate(mesh)

# a measure of cell size
h = sqrt(CellVolume(mesh))

# setting up vertical direction
y_abs = sqrt(y**2)
yhat = as_vector((0, y)) / y_abs

# Global Constants:
max_num_timesteps = 200
simu_time = 0.0

# Stokes related constants:
Ra = Constant(1e6)  # Rayleigh Number

# Temperature related constants:
delta_t = Constant(2e-6)  # Time-step
kappa = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson
# scheme where theta_ts = 0.5.
theta_ts = 0.5


# Print function to ensure log output is only written
# on processor zero (if running in parallel) ####
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

# Geometry and Spatial Discretization:
# Set up function spaces - currently using the P2P1 element pair :
V = VectorFunctionSpace(mesh, "CG", 2)  # Velocity function space (vector)
W = FunctionSpace(mesh, "CG", 1)  # Pressure function space (scalar)
Q = FunctionSpace(mesh, "CG", 2)  # Temperature function space (scalar)

# Set up mixed function space and associated test functions:
Z = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)
Y = TestFunction(Q)

# Set up fields on these function spaces - split into each 
# component so that they are easily accessible:
z = Function(Z)  # a field over the mixed function space Z.
u, p = split(z)     # can we nicely name mixed function space fields?

# reference velocity function
# to be used in the inversion
ref_vel = Function(V, name="Reference_Velocity")

target_cfl_no = 2.5
max_timestep = 1.00
maximum_timestep = 0.1
increase_tolerance = 1.5
simu_time = 0.0


# Timestepping - CFL related stuff:
def compute_timestep(u, current_delta_t):
    """Return the timestep, based upon the CFL criterion"""

    ref_vel.interpolate(dot(JacobianInverse(mesh), u))
    ts_min = 1. / mesh.comm.allreduce(ref_vel.dat.data.max(), MPI.MAX)
    # Grab (smallest) maximum permitted on all cores:
    ts_max = min(float(current_delta_t) * increase_tolerance, maximum_timestep)
    # Compute timestep:
    tstep = min(ts_min * target_cfl_no, ts_max)
    log(f"ts_max= {ts_max}")
    return tstep


# Temperature
T_old = Function(Q, name="Temperature")
T_old.interpolate(0.5*(erf((1-X[1])*thickness_val)+erf(-X[1]*thickness_val)+1) + 0.1*exp(-0.5*((X-as_vector((0.5, 0.2)))/Constant(0.1))**2))

# Defining temperature field and initialise it with old temperature
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
F_stokes += - (dot(N, yhat) * Ra * T_theta) * dx
F_stokes += - div(u) * M * dx

# Setting free-slip BC for top and bottom
bcu_topbase = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id))
bcu_rightleft = DirichletBC(Z.sub(0).sub(0), 0.0, (left_id, right_id))
all_bcu = [bcu_topbase, bcu_rightleft]

# Setting Dirihlet BC for top and bottom of the temperature field
bct_top = DirichletBC(Q, 0.0, (top_id))
bct_base = DirichletBC(Q, 1.0, (bottom_id))

# Pressure nullspace
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0),
                                    VectorSpaceBasis(constant=True)]
                                    )

# Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx\
            + Y * dot(u, grad(T_theta)) * dx\
            + dot(grad(Y), kappa * grad(T_theta)) * dx


# split of the fields to access for IO
u_, p_ = z.split()
u_.rename('Velocity')
p_.rename('Pressure')

# Printing out the degrees of freedomG
log('global number of nodes P1 coeffs/nodes:', W.dim())

# Setup problem and solver objects so we can reuse (cache) solver setup
stokes_problem = NonlinearVariationalProblem(F_stokes, z, bcs=all_bcu)
stokes_solver = NonlinearVariationalSolver(stokes_problem,  solver_parameters=solver_parameters,
     nullspace=p_nullspace, transpose_nullspace=p_nullspace)

energy_problem = NonlinearVariationalProblem(F_energy, T_new, bcs=[bct_base, bct_top])
energy_solver = NonlinearVariationalSolver(energy_problem, solver_parameters=solver_parameters)


# Write functions out in VTK format
state_vtu_file = File('visual/state.pvd')

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems,
    # where everything is on the LHS (as above)
    # and the RHS == 0.
    stokes_solver.solve()

    # Write output:
    if timestep % 10 == 0:
        log(f"Output: {simu_time:.3e}, {timestep:.0f}")
        state_vtu_file.write(u_, p_, T_new)

    new_delta_t = compute_timestep(u, delta_t)

    # Temperature system:
    energy_solver.solve()

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log((f"idx: {timestep} "
         f"time: {simu_time:.1e} "
         f"del_t: {delta_t.__float__():.1e}"
         )
        )


# Output the last state
state_vtu_file.write(u_, p_, T_new)

# Generating the reference temperature field for the adjoint
checkpoint_data = CheckpointFile("Final_State.h5", "w")
checkpoint_data.save_mesh(mesh)
checkpoint_data.save_function(T_new)
checkpoint_data.close()
