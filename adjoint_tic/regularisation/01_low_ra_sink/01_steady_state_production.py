"""
    Generation of forward problem - for an adjoint simulation
    The test problem is no heat-flow at boundaries, with free-slip on all surfaces
"""
from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

# Quadrature degree: 
dx = dx(degree=6) 

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
max_num_timesteps      = 20000
target_cfl_no          = 2.5
max_timestep           = 1.00
maximum_timestep    = 0.1
increase_tolerance = 1.5
simu_time = 0.0
# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number

# Temperature related constants:
delta_t                = Constant(1e-4) # Time-step
kappa                  = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson scheme where theta_ts = 0.5:
theta_ts               = 0.5

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
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
         "ksp_rtol": 1e-4,
         #"ksp_converged_reason": None,
         "pc_type": "python",
         "pc_python_type": "firedrake.AssembledPC",
         "assembled_pc_type": "gamg",
         "assembled_pc_gamg_threshold": 0.01,
         "assembled_pc_gamg_square_graph": 100,
     },
    "fieldsplit_1": {
        "ksp_type": "fgmres",
        "ksp_rtol": 1e-3,
        #"ksp_converged_reason": None,
        "pc_type": "python",
        "pc_python_type": "firedrake.MassInvPC",
        "Mp_ksp_rtol": 1e-5,
        "Mp_ksp_type": "cg",
        "Mp_pc_type": "sor",
    }
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

# Timestepping - CFL related stuff:
ref_vel = Function(V, name="Reference_Velocity")

def compute_timestep(u, current_delta_t):
    """Return the timestep, based upon the CFL criterion"""

    ref_vel.interpolate(dot(JacobianInverse(mesh), u))
    ts_min = 1. / mesh.comm.allreduce(ref_vel.dat.data.max(), MPI.MAX)
    # Grab (smallest) maximum permitted on all cores:
    ts_max = min(float(current_delta_t) * increase_tolerance, maximum_timestep)
    # Compute timestep:
    tstep = min(ts_min * target_cfl_no, ts_max)
    return tstep

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Set up temperature field and initialise based upon coordinates:
T_old    = Function(Q, name="OldTemperature")

# Having a single hot blob as initial condition
blb_ctr_h = as_vector((0.0, 0.00)) 
blb_gaus = Constant(0.1)

# A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
T_old.interpolate(1.0 - y + 0.1*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2));

# Defining temperature field and initialise it with old temperature
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

# Setting free-slip BC for top and bottom
bcu_topbase     = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id))
bcu_rightleft   = DirichletBC(Z.sub(0).sub(0), 0.0, (left_id, right_id))


# Setting Dirihlet BC for top and bottom of the temperature field
bct_top  = DirichletBC(Q, 0.0, (top_id))
bct_base = DirichletBC(Q, 1.0, (bottom_id))

# Pressure nullspace                                                                                                                                                                                                                                                                      
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

# Generating near_nullspaces for GAMG:
z_rotV = Function(V).interpolate(as_vector((-X[1], X[0])))
nns_x = Function(V).interpolate(Constant([1., 0.]))
nns_y = Function(V).interpolate(Constant([0., 1.]))
V_near_nullspace = VectorSpaceBasis([nns_x, nns_y, z_rotV])
V_near_nullspace.orthonormalize()
Z_near_nullspace = MixedVectorSpaceBasis(Z, [V_near_nullspace, Z.sub(1)])

def absv(u):
    """Component-wise absolute value of vector for SU stabilisation"""
    return as_vector([abs(ui) for ui in u])

def beta(Pe):
    """Component-wise beta formula Donea and Huerta (2.47a) for SU stabilisation"""
    return as_vector([1/tanh(Pei+1e-6) - 1/(Pei+1e-6) for Pei in Pe])

# copying SU scheme from global gplates
J = Function(TensorFunctionSpace(mesh, 'DQ', 1), name='Jacobian').interpolate(Jacobian(mesh))
Pe = absv(dot(u, J)) / 2
nubar = dot(Pe, beta(Pe))
Y_SU = Y + nubar / dot(u, u) * dot(u, grad(Y))

### Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y_SU*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# Write output in VTK format:
state_vtu_file = File('vis_steadystate/state.pvd')

# For some reason this only works here!!!
u_, p_    = z.split() 
u_.rename('Velocity') 
p_.rename('Pressure')

# Printing out the degrees of freedomG
log('global number of nodes P1 coeffs/nodes:', W.dim())


# Setup problem and solver objects so we can reuse (cache) solver setup
stokes_problem = NonlinearVariationalProblem(F_stokes, z, bcs=[bcu_topbase, bcu_rightleft])
stokes_solver = NonlinearVariationalSolver(stokes_problem, solver_parameters=None,
     appctx={"mu": mu}, nullspace=p_nullspace, transpose_nullspace=p_nullspace,
     near_nullspace=Z_near_nullspace
)

energy_problem = NonlinearVariationalProblem(F_energy, T_new, bcs=[bct_base, bct_top])
energy_solver = NonlinearVariationalSolver(energy_problem, solver_parameters=Noneenergy_iterative)

## Checkpointing 
#u_checkpoint = CheckpointFile("velocities.h5", 'w')
#u_checkpoint.save_mesh(mesh)

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    stokes_solver.solve()

    # writing out the velocity field, which will be used for reconstructions
    #u_checkpoint.save_function(u, idx=timestep)

    delta_t.assign(compute_timestep(u, delta_t))

    # Write output:
    if timestep % 200 == 0:
        log("Output:", simu_time, timestep)
        state_vtu_file.write(u_, p_, T_new)
        # Generating the reference temperature field for the adjoint
        #checkpoint_data = CheckpointFile(f"State_{timestep}", "w")
        #checkpoint_data.save_mesh(mesh)
        #checkpoint_data.save(T_new)
        #checkpoint_data.close() 

    # Temperature system:
    energy_solver.solve()

    # updating time
    simu_time += float(delta_t)

    # computing how much temperature has changed
    T_change = sqrt(assemble((T_old-T_new)**2 * dx))

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log(f"Timestep Number: {timestep} Timestep: {delta_t.__float__():.2e} T_change {T_change:.3e}")

    if T_change<=1e-5:
        break

#u_checkpoint.close()

checkpoint_data = CheckpointFile(f"State_SteadyState.h5", "w")
checkpoint_data.save_mesh(mesh)
checkpoint_data.save_function(T_new)
checkpoint_data.close()
