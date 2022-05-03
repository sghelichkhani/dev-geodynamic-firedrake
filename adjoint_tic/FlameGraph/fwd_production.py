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
steady_state_tolerance  = 1e-7
max_num_timesteps       = 100
target_cfl_no           = 1.0
maximum_timestep        = 1.00
increase_tolerance	= 1.5
time                    = 0.0


# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number

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

# Timestepping - we need ref_vel to map u
ref_vel = Function(V, name="Reference_Velocity")

def compute_timestep(u, current_delta_t): 
    """Return the timestep, based upon the CFL criterion"""
    ref_vel.interpolate(dot(JacobianInverse(mesh), u))
    ts_min = 1. / mesh.comm.allreduce(ref_vel.dat.data.max(), MPI.MAX)
    # Grab (smallest) maximum permitted on all cores:
    ts_max = min(float(current_delta_t)*increase_tolerance, maximum_timestep)
    # Compute timestep:
    tstep = min(ts_min*target_cfl_no, ts_max)
    return tstep

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Set up temperature field and initialise from provided file 
T_old    = Function(Q, name="OldTemperature")

# Loading initial condition from previously generated files 
InitialStateFile = DumbCheckpoint("./InitialCondition", mode=FILE_READ)
InitialStateFile.load(T_old, "Temperature")
InitialStateFile.close()

# Defining temperature field and initialise it with old temperature
T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

# helper function to compute horizontal layer averages
Tlayer = Function(Qlayer, name='LayerTemp')  # stores values of temp in one layer
Tavg = Function(Q, name='LayerAveragedTemp')  # averaged temp function returned by function
Rmin_area = assemble(Constant(1.0, domain=mesh1d)*dx)  # area of CMB

def layer_average(T):
    vnodes = disc_n*2 + 1  # n/o Q2 nodes in the vertical
    hnodes = Qlayer.dim()  # n/o Q2 nodes in each horizontal layer
    assert hnodes*vnodes == Q.dim()
    for i in range(vnodes):
        Tlayer.dat.data[:] = T.dat.data_ro[i::vnodes]
        # NOTE: this integral is performed on mesh2d, which always has r=Rmin, but we normalize
        Tavg.dat.data[i::vnodes] = assemble(Tlayer*dx) / Rmin_area
    return Tavg

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
bcu_topbase     = DirichletBC(Z.sub(0), 0.0, (top_id, bottom_id))
bcu_rightleft   = DirichletBC(Z.sub(0), 0.0, (left_id, right_id))

# Setting Dirichlet temperature boundary condition
bcT_top     = DirichletBC(Q, 0.0, (top_id))
bcT_base    = DirichletBC(Q, 1.0, (bottom_id))

# Pressure nullspace
p_nullspace = MixedVectorSpaceBasis(Z, [Z.sub(0), VectorSpaceBasis(constant=True)])

### Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# Write output files in VTK format:
u_file = File('Visual/velocity.pvd')
p_file = File('Visual/pressure.pvd')
t_file = File('Visual/temperature.pvd')

# For some reason this only works here!!!
u, p    = z.split() 
u.rename('Velocity') 
p.rename('Pressure')

# Printing out the degrees of freedom 
log('global number of nodes P1 coeffs/nodes:', W.dim())

# A simulation time to track how far we are
simu_time = 0.0

# Stokes Solver
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

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    stokes_solver.solve()

    ## Current time-stepping
    #current_delta_t = delta_t

    #delta_t.assign(compute_timestep(u, current_delta_t))  # Compute adaptive time-step
    time += float(delta_t)

    # Write output:
    if timestep % 5 == 0:
        Tavg = layer_average(T_new)
        log("Output:", simu_time, timestep)
        u_file.write(u)
        p_file.write(p)
        t_file.write(T_new, Tavg)

        # Generating reference states 
        checkpoint_data = DumbCheckpoint(str("State_%2.2i" %(timestep)),\
                                         single_file=True, mode=FILE_CREATE, comm=mesh.comm)
        checkpoint_data.store(T_new)
        checkpoint_data.store(Tavg)
        checkpoint_data.close()

    # Temperature system:
    energy_solver.solve()

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log(f"Timestep: {timestep}, {float(delta_t)}")

# Because here we are one step further
timestep += 1

# Generating reference states 
checkpoint_data = DumbCheckpoint(str("State_%2.2i" %(timestep)),\
                                 single_file=True, mode=FILE_CREATE, comm=mesh.comm)
checkpoint_data.store(T_new)
checkpoint_data.store(Tavg)
checkpoint_data.close()

