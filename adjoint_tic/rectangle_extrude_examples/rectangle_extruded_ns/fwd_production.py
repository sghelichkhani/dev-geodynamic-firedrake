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

# setting up spatial coordinates
X  =  x, y = SpatialCoordinate(mesh)

# a measure of cell size
h	  = sqrt(CellVolume(mesh))

# setting up vertical direction
y_abs     = sqrt(y**2)
yhat  = as_vector((0,y)) / y_abs

# Global Constants:
steady_state_tolerance = 1e-7
max_num_timesteps      = 65
target_cfl_no          = 2.5
max_timestep           = 1.00

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number

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

# Timestepping - CFL related stuff:
ts_func = Function(Q) # Note that time stepping should be dictated by Temperature related mesh.

def compute_timestep(u):
    # A function to compute the timestep, based upon the CFL criterion
    ts_func.interpolate( h / sqrt(dot(u,u)))
    ts_min = ts_func.dat.data.min()
    ts_min = mesh.comm.allreduce(ts_min, MPI.MIN)
    #return min(ts_min*target_cfl_no,max_timestep)
    return 5.0e-05 

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Set up temperature field and initialise based upon coordinates:
T_old    = Function(Q, name="OldTemperature")

# Having a single hot blob on 1.5, 0.0
blb_ctr_h = as_vector((0.5, 0.85)) 
blb_gaus = Constant(0.04)

# A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
T_old.interpolate(0.5 - 0.3*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2));

# Defining temperature field and initialise it with old temperature
T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

#********************** Set-up Equations ******************** 

### Initially deal with Stokes equations ###
# Equation in weak (ufl) form 
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

### Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

## Prescribed temperature for top and bottom
#bct_base = DirichletBC(Q, 1.0, bottom_id)
#bct_top  = DirichletBC(Q, 0.0, top_id)

# Write output files in VTK format:
u_file = File('FWDREFmodel/velocity.pvd')
p_file = File('FWDREFmodel/pressure.pvd')
t_file = File('FWDREFmodel/temperature.pvd')

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
stokes_solver  = LinearVariationalSolver(stokes_problem, solver_parameters=solver_parameters)

q_tri = TrialFunction(Q)
F_energy_lin = replace(F_energy, {T_new:q_tri})
a_energy, L_energy = lhs(F_energy_lin), rhs(F_energy_lin)
energy_problem = LinearVariationalProblem(a_energy, L_energy, T_new)
energy_solver  = LinearVariationalSolver(energy_problem, solver_parameters=solver_parameters)

# Now perform the time loop:
for timestep in range(0, max_num_timesteps):
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    stokes_solver.solve() 

    # updating time-step based on velocities
    delta_t.assign(compute_timestep(u)) # Compute adaptive time-step

    # Temperature system:
    energy_solver.solve()

    # updating time
    simu_time += float(delta_t)

    # Write output:
    if timestep % 5 == 0:
        log("Output:", simu_time, timestep)
        u_file.write(u)
        p_file.write(p)
        t_file.write(T_new)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log("Timestep Number: ", timestep, " Timestep: ", float(delta_t))

# Generating the reference temperature field for the adjoint
checkpoint_data = DumbCheckpoint("final_state", single_file=True, mode=FILE_CREATE, comm=mesh.comm)
checkpoint_data.store(T_new)
checkpoint_data.close()

log("Final time:", simu_time)
u_file.write(u)
p_file.write(p)
t_file.write(T_new)
