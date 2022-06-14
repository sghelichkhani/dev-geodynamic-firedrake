""""
    Generation of forward problem - for an adjoint simulation
    The test problem is no heat-flow at boundaries, with free-slip on all surfaces
"""
import matplotlib.pyplot as plt
from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc
from lib_averaging import LayerAveraging

# defining a mesh to define average
x_max = 1.0
#  how many intervals along x/y directions 
disc_n = 150
mesh1d = IntervalMesh(disc_n, length_or_left=0.0, right=x_max) 

# Some important constants etc
# Quadrature degree:
dx = dx(degree=6)

# logging.set_log_level(1)
# logging.set_level(1)

# Loading the initial State
with CheckpointFile("Initial_State.h5", "r") as t_ic_checkpoint:
    mesh = t_ic_checkpoint.load_mesh("firedrake_default_extruded")
    T_old = t_ic_checkpoint.load_function(mesh, "Temperature")


# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'
left_id, right_id = 1, 2

# spatial coordinates
X = x, y = SpatialCoordinate(mesh)

# a measure of cell size
h = sqrt(CellVolume(mesh))

# setting up vertical direction
y_abs = sqrt(y**2)
yhat = as_vector((0, y)) / y_abs

# Global Constants:
max_num_timesteps = 150
simu_time = 0.0

# Stokes related constants:
Ra = Constant(1e6)  # Rayleigh Number

# Temperature related constants:
delta_t = Constant(1e-6)  # Time-step
kappa = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson
# scheme where theta_ts = 0.5.
theta_ts = 0.5


# Print function to ensure log output is only written
# on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args)

# Energy Equation Solver Parameters:
energy_iterative = {
    "mat_type": "aij",
    "snes_type": "ksponly",
    "ksp_type": "gmres",
    "ksp_rtol": 1e-5,
    #"ksp_converged_reason": None,
    "pc_type": "sor",
}

stokes_iterative = {
     "mat_type": "matfree",
     "snes_type": "ksponly",
     "ksp_type": "preonly",
     #"ksp_converged_reason": None,
     "pc_type": "fieldsplit",
     "pc_fieldsplit_type": "schur",
     "pc_fieldsplit_schur_type": "full",
     "fieldsplit_0": {
         "ksp_type": "cg",
         "ksp_rtol": 1e-5,
         #"ksp_converged_reason": None,
         "pc_type": "python",
         "pc_python_type": "firedrake.AssembledPC",
         "assembled_pc_type": "gamg",
         "assembled_pc_gamg_threshold": 0.01,
         "assembled_pc_gamg_square_graph": 100,
     },
    "fieldsplit_1": {
        "ksp_type": "fgmres",
        "ksp_rtol": 1e-4,
        #"ksp_converged_reason": None,
        "pc_type": "python",
        "pc_python_type": "firedrake.MassInvPC",
        "Mp_ksp_rtol": 1e-5,
        "Mp_ksp_type": "cg",
        "Mp_pc_type": "sor",
    }
}


# Geometry and Spatial Discretization:
# Set up function spaces - currently using the P2P1 element pair :
V = VectorFunctionSpace(mesh, "CG", 2)  # Velocity function space (vector)
W = FunctionSpace(mesh, "CG", 1)  # Pressure function space (scalar)
Q = FunctionSpace(mesh, "CG", 2)  # Temperature function space (scalar)
Qlayer = FunctionSpace(mesh1d, "CG", 2)  # Temperature function space on the 1D mesh 

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

# Timestepping - CFL related stuff:
def compute_timestep(u, current_delta_t):
    """Return the timestep, based upon the CFL criterion"""

    ref_vel.interpolate(dot(JacobianInverse(mesh), u))
    ts_min = 1. / mesh.comm.allreduce(ref_vel.dat.data.max(), MPI.MAX)
    # Grab (smallest) maximum permitted on all cores:
    ts_max = min(float(current_delta_t) * increase_tolerance, maximum_timestep)
    # Compute timestep:
    tstep = min(ts_min * target_cfl_no, ts_max)
    return tstep


# helper function to compute horizontal layer averages
Tlayer = Function(Qlayer, name='LayerTemp')  # stores values of temp in one layer

# T advection diffusion equation Prerequisites:
T_average = Function(Q, name="OneDimTemperature")


# A layer average to define temperature
def layer_average(T):
    vnodes = disc_n*2 + 1  # n/o Q2 nodes in the vertical
    hnodes = Qlayer.dim()  # n/o Q2 nodes in each horizontal layer
    assert hnodes*vnodes == Q.dim()
    for i in range(vnodes):
        Tlayer.dat.data[:] = T.dat.data_ro[i::vnodes]
        # NOTE: this integral is performed on mesh2d
        T_average.dat.data[i::vnodes] = assemble(Tlayer*dx) 
    return T_average

# Temperature
T_deviatoric = Function(Q, name="TemperatureDev")

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
# Generating near_nullspaces for GAMG:
z_rotV = Function(V).interpolate(as_vector((-X[1], X[0])))
nns_x = Function(V).interpolate(Constant([1., 0.]))
nns_y = Function(V).interpolate(Constant([0., 1.]))
V_near_nullspace = VectorSpaceBasis([nns_x, nns_y, z_rotV])
V_near_nullspace.orthonormalize()
Z_near_nullspace = MixedVectorSpaceBasis(Z, [V_near_nullspace, Z.sub(1)])

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
stokes_solver = NonlinearVariationalSolver(stokes_problem, solver_parameters=stokes_iterative,
     appctx={"mu": mu}, nullspace=p_nullspace, transpose_nullspace=p_nullspace,
     near_nullspace=Z_near_nullspace
)

energy_problem = NonlinearVariationalProblem(F_energy, T_new, bcs=[bct_base, bct_top])
energy_solver = NonlinearVariationalSolver(energy_problem, solver_parameters=energy_iterative)

# Checkpointing
u_tave_checkpoint = CheckpointFile("Ref_velocities.h5", 'w')
u_tave_checkpoint.save_mesh(mesh)


# Compute the average radial temperature profile
layer_average(T_new)

# Write the average temperature out
# This we can use for regularisation
u_tave_checkpoint.save_function(T_average, idx=0)

# Write functions out in VTK format
state_vtu_file = File('vis_ref_simulation/state.pvd')

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems,
    # where everything is on the LHS (as above)
    # and the RHS == 0.
    stokes_solver.solve()

    # writing out the velocity field, which will be used for reconstructions
    u_tave_checkpoint.save_function(u_, idx=timestep)

    # compute average
    layer_average(T_new)

    # Write output:
    if timestep % 10 == 0:
        log(f"Output: {simu_time:.3e}, {timestep:.0f}")
        T_deviatoric.interpolate(T_new - T_average)
        state_vtu_file.write(u_, p_, T_deviatoric)

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

u_tave_checkpoint.close()

# compute the average temperature
layer_average(T_new)

# Output the last state
T_deviatoric.interpolate(T_new - T_average)
state_vtu_file.write(u_, p_, T_deviatoric)

# Generating the reference temperature field for the adjoint
checkpoint_data = CheckpointFile("Final_State.h5", "w")
checkpoint_data.save_mesh(mesh)
checkpoint_data.save_function(T_new)
checkpoint_data.save_function(T_average)
checkpoint_data.close()
