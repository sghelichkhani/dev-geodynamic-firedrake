from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

x_max = 2.0; y_max = 1.0;
x_dis = 160; y_dis = 80;


mesh2d = IntervalMesh(ncells=x_dis, length_or_left=0., right=x_max)
mesh = ExtrudedMesh(mesh2d, layers=y_dis, layer_height=y_max/y_dis)


# Mesh and associated physical boundary IDs:
#mesh    = PeriodicRectangleMesh(x_dis, y_dis, x_max, y_max, direction='x')
bottom_id, top_id = 'bottom', 'top'  # with PeriodicRectangleMesh: 1,2 for bottom and top
left_id, right_id = 1, 2  # with PeriodicRectangleMesh: 1,2 for bottom and top

# Global Constants:
steady_state_tolerance = 1e-8
max_timesteps          = 500 
model_time             = 0.0
max_timestep           = 1.0

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number
k                      = Constant((0,1))
X = x,y                = SpatialCoordinate(mesh)
h                      = sqrt(CellVolume(mesh))

# Temperature related constants:
delta_t                = Constant(1.0e-9) # Initial time-step
kappa                  = Constant(1.0)  # Thermal diffusivity
target_cfl_no          = 2.5
# Temporal discretisation - Using a Crank-Nicholson scheme where theta_ts = 0.5:
theta_ts               = 0.5

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args)

# Stokes Equation Solver Parameters:
stokes_solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_type': 'mumps',
    'mat_type': 'aij'
}

# Temperature Equation Solver Parameters:
temperature_solver_parameters = {
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
Z    = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)
Y    = TestFunction(Q)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z    = Function(Z) # a field over the mixed function space Z.
u, p = split(z) # can we nicely name mixed function space fields?

# time step function
ts_func = Function(Q)

def compute_timestep(u):
    # A function to compute the timestep, based upon the CFL criterion
    ts_func.interpolate( h / sqrt(dot(u,u)))
    ts_min = ts_func.dat.data.min()
    ts_min = mesh.comm.allreduce(ts_min, MPI.MIN)
    return min(ts_min*target_cfl_no,max_timestep)


#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# The gaussian structure
blb_ctr_h = as_vector((1.5, 0.95))
blb_gaus = Constant(0.1)
# 

# Set up temperature field and initialise based upon coordinates:
T_old   = Function(Q, name="OldTemperature")
T_old.project((y_max - y) - 0.3*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2))


T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

#########################################################################################################
############################################ Setup Equations ############################################
#########################################################################################################

### Initially deal with Stokes equations ###

# Equation in weak (ufl) form - note that continuity equation is added here - need to better understand why:
# Set up in residual form to ensure non-linear solvers are used.
mu_field    = Function(W, name="Viscosity")
mu          = 1.0*exp(-ln(10)*T_new) # Variable viscosity
tau         = mu * (grad(u)+transpose(grad(u))) # Strain-rate tensor:
tauii       = sqrt((tau[0,0]**2 + tau[0,1]**2 + tau[0,1]**2 + tau[1,1]**2) / 2.)
F_stokes    = inner(grad(N), tau) * dx + dot(N,grad(p)) * dx - (dot(N,k)*Ra*T_theta) * dx
F_stokes   += dot(grad(M),u) * dx # Continuity equation

# Set up boundary conditions for Stokes: We first need to extract the velocity field from the mixed
# function space Z. This is done using Z.sub(0). We subsequently need to extract the x and y components
# of velocity. This is done using an additional .sub(0) and .sub(1), respectively. Note that the final arguments
# here are the physical boundary ID's.
#bcv_val_top  = Function(V).interpolate(as_vector((2000.0, 0.0)))
#bcv_top      = DirichletBC(Z.sub(0), bcv_val_top, top_id) 
#bcv_top      = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id))
#bcv_val_base = Function(V).interpolate(as_vector((0.0, 0.0)))
#bcv_base     = DirichletBC(Z.sub(0), bcv_val_base, bottom_id)
bcv_fs_tb      = DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id))
bcv_fs_lr      = DirichletBC(Z.sub(0).sub(0), 0.0, (left_id, right_id))

### Next deal with Temperature advection-diffusion equation: ###
delta_x      = sqrt(CellVolume(mesh))
Y_SUPG       = Y + dot(u,grad(Y)) * (delta_x / (2*sqrt(dot(u,u))))
F_energy     = Y_SUPG * ((T_new - T_old) / delta_t) * dx + Y_SUPG*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx
bct_base     = DirichletBC(Q, 1.0, bottom_id)
bct_top      = DirichletBC(Q, 0.0, top_id)


# Write output files in VTK format:
u, p = z.split() # Do this first to extract individual velocity and pressure fields:
u._name = 'Velocity' # This doesn't seem to work!
p._name = 'Pressure'
u_file     = File('velocity.pvd')
p_file     = File('pressure.pvd')
t_file     = File('temperature.pvd')
mu_file    = File('viscosity.pvd')


# Starting O file for storing Temperatures
checkpoint_evolution = DumbCheckpoint("./ThermalCheckpoints",mode=FILE_CREATE)

# Now perform the time loop:
for timestep in range(0,max_timesteps):
    log(f"Timestep Number: {timestep}, Time: {model_time}, delta_t: {float(delta_t)}")

    if(timestep % 10 == 0):
        # Write output:
        u_file.write(u)
        p_file.write(p)
        t_file.write(T_new)
        mu_field.interpolate(mu)
        mu_file.write(mu_field)

        # storing the evolution
        checkpoint_evolution.set_timestep(t=model_time, idx=timestep)
        checkpoint_evolution.store(T_new)        

    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0.
    solve(F_stokes==0, z, bcs=[bcv_fs_tb, bcv_fs_lr], solver_parameters=stokes_solver_parameters)
    
    delta_t.assign(compute_timestep(u))

    # Temperature system:
    solve(F_energy==0, T_new, bcs=[bct_base,bct_top], solver_parameters=temperature_solver_parameters)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Update model_time
    model_time += float(delta_t)

checkpoint_evolution.close()

