from firedrake import *
from firedrake.petsc import PETSc
from mpi4py import MPI

# Set up geometry:
rmin, rmax      = 1.22, 2.22
ncells, nlayers = 256, 64

# Construct a circle mesh and then extrude into a cylinder:
mesh1d = CircleManifoldMesh(ncells, radius=rmin)
mesh   = ExtrudedMesh(mesh1d, layers=nlayers, extrusion_type='radial')
# Physical boundary IDs from extruded mesh:
bottom_id, top_id = 'bottom', 'top'

## Constants and set-up P2 version (popped) of mesh
P1    = FunctionSpace(mesh, "CG", 1)
V     = VectorFunctionSpace(mesh, "CG", 2)
x, y  = SpatialCoordinate(mesh)
r     = sqrt(x**2 + y**2)
r_p1  = interpolate(r, P1)
xy_p2 = interpolate(as_vector((x,y))/r*r_p1, V)
mesh  = Mesh(xy_p2)

# now redefine everything based on the super/iso parametric P2 mesh
x, y  = SpatialCoordinate(mesh)
n     = FacetNormal(mesh)
h     = sqrt(CellVolume(mesh))
r     = sqrt(x**2 + y**2)
k     = as_vector((x, y)) / r
domain_volume = assemble(1.*dx(domain=mesh))

# Parameters relating to time-stepping:
steady_state_tolerance = 1e-8
max_timesteps          = 25000
target_cfl_no          = 2.5
maximum_timestep       = 1.0
increase_tolerance     = 1.2
time                   = 0.0
dim                    = 2

# Stokes related constants (note that since these are included in UFL, they are wrapped inside Constant):
mu      = Constant(1.0)   # Viscosity - constant for this isoviscous case.

Ra      = Constant(1e4)   # Rayleigh number.
C_ip    = Constant(20)    # Fudge factor for interior penalty term used in weak imposition of BCs.

# Temperature equation related constants:
delta_t = Constant(1e-6) # Initial time-step.
kappa   = Constant(1.0)  # Thermal diffusivity.

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    PETSc.Sys.Print(*args)

#### File for logging output diagnostics through simulation.
def log_params(f, str):
    if mesh.comm.rank == 0:
        f.write(str + "\n")
        f.flush()

#########################################################################################################
######################################### Solver parameters:  ###########################################
#########################################################################################################

# Stokes Equation Solver Parameters:
stokes_solver_parameters = {
    "mat_type": "matfree",
    "snes_type": "ksponly",
    "ksp_type": "preonly",
    "pc_type": "fieldsplit",
    "pc_fieldsplit_type": "schur",
    "pc_fieldsplit_schur_type": "full",
    "fieldsplit_0": {
        "ksp_type": "cg",
        "pc_type": "python",
        "pc_python_type": "firedrake.AssembledPC",
        "assembled_pc_type": "gamg",
        "assembled_pc_gamg_threshold_scale": 1.0,
        "assembled_pc_gamg_threshold": 0.01,
        "assembled_pc_gamg_coarse_eq_limit": 800,
        # Note the following 5 options are required only when there is a velocity null space (which is the case here)
        "assembled_mg_coarse_ksp_type": "preonly",
        "assembled_mg_coarse_pc_type": "sor",
        "assembled_mg_coarse_ksp_max_it": 10,
        "assembled_mg_coarse_ksp_rtol": 0., # we always want 10 iterations
        "assembled_mg_coarse_ksp_atol": 0., # we always want 10 iterations
        "ksp_rtol": 1e-7,     
        "ksp_converged_reason": None,
    },
    "fieldsplit_1": {
        "ksp_type": "fgmres",
        "ksp_converged_reason": None,
        "pc_type": "python",
        "pc_python_type": "firedrake.MassInvPC",
        "Mp_ksp_type": "cg",
        "Mp_pc_type": "sor",
        "ksp_rtol": 1e-6,
    }
}

# Temperature Equation Solver Parameters:
temperature_solver_parameters = {
        "snes_type": "ksponly",
        "ksp_rtol": 1e-5,
        "ksp_type": "gmres",
        "pc_type": "sor",
        "mat_type": "aij"
}

#########################################################################################################
######################################## Spatial Discretization: ########################################
#########################################################################################################

# Set up function spaces - currently using the bilinear Q2Q1 element pair:
V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W    = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Q    = FunctionSpace(mesh, "CG", 2) # Temperature function space (scalar)

log("Number of Velocity DOF:", V.dim()*3)
log("Number of Pressure DOF:", W.dim())
log("Number of Velocity and Pressure DOF:", V.dim()*3+W.dim())
log("Number of Temperature DOF:", Q.dim())

# Set up mixed function space and associated test functions:
Z    = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)
Y    = TestFunction(Q)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z    = Function(Z) # a field over the mixed function space Z.
u, p = split(z) # returns a symbolic UFL expression for velocity (u) and pressure (p)

#########################################################################################################
############################################ Time Stepping: #############################################
#########################################################################################################

# Timestepping - CFL related stuff:
delta_x = sqrt(CellVolume(mesh))
ts_func = Function(Q) # Note that time stepping should be dictated by Temperature related mesh.

def compute_timestep(u, current_delta_t):
    # A function to compute the timestep, based upon the CFL criterion.
    ts_func.interpolate(delta_x / sqrt(dot(u,u)))
    ts_min = ts_func.dat.data.min()
    ts_min = mesh.comm.allreduce(ts_min, MPI.MIN)
    ts_max = min(current_delta_t.dat.data[0]*increase_tolerance, maximum_timestep)
    return max(ts_min*target_cfl_no, ts_max)

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Set up temperature field and initialise:
T_old   = Function(Q, name="OldTemperature")
T_old.interpolate(rmax-r + 0.02*cos(4.*atan_2(y,x))*sin((r-rmin)*pi))

# To start set T_new to T_old.
T_new   = Function(Q, name="Temperature")
T_new.assign(T_old)

# Temporal discretisation - Using a Crank-Nicholson scheme where theta = 0.5:
theta   = 1.0
T_theta = theta * T_new + (1-theta) * T_old

#########################################################################################################
############################################ Setup Equations ############################################
#########################################################################################################

### Initially deal with Stokes equations ###
#mu  = exp(-ln(1e4)*T_theta)

# deviatoric stresses
def tau(u): return  mu * (grad(u)+transpose(grad(u)))

# traction field 
def trac(u,p): return dot(tau(u),n) - p*n

# nitsche free slip BCs
# (p+1)*(p+d)/d; with p:= polynomial degree, d:= dimension => (2+1)(2+2)/2
nitsche_fs  = - dot(N,n)*dot(n,trac(u,p))*ds_tb - dot(u,n)*dot(n,trac(N,M))*ds_tb + C_ip*((2+1)*(2+dim)/dim)*FacetArea(mesh)/CellVolume(mesh)*dot(u,n)*dot(N,n)*ds_tb

# UFL for Stokes equations (note that given we are imposing free-slip BCs weakly, we have done some integration by parts:
F_stokes  = inner(grad(N), tau(u)) * dx -div(N) * p * dx - (dot(N,k)*Ra*T_theta) * dx
F_stokes += -div(u) * M * dx # Continuity equation
F_stokes += nitsche_fs

# Nullspaces and near-nullspaces:
rotZ     = Function(Z)
rotV, _  = rotZ.split()
rotV.interpolate(as_vector((-y, x)))
rotVSpace = VectorSpaceBasis([rotV])
rotVSpace.orthonormalize()

# Constant nullspace for pressure
p_nullspace = VectorSpaceBasis(constant=True)

# Setting the mixed nullspace
Z_nullspace = MixedVectorSpaceBasis(Z, [rotVSpace, p_nullspace])

# Generating near_nullspaces for GAMG, starting with constant modes:
c0V = Function(V)
c0V.interpolate(Constant([1., 0.]))
c1V = Function(V)
c1V.interpolate(Constant([0., 1.]))
# Rotation around Z axis (using rotV from above here messes with nullspaces due to orthonormalisation)
cr0V = Function(V)
cr0V.interpolate(as_vector((-y, x)))

V_near_nullspace = VectorSpaceBasis([cr0V, c0V, c1V])
V_near_nullspace.orthonormalize()
Z_near_nullspace = MixedVectorSpaceBasis(Z, [V_near_nullspace, Z.sub(1)])

# Next deal with Temperature advection-diffusion equation - in this example, we do not use stabilisation:
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y * dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# Temperature boundary conditions
bct_base = DirichletBC(Q, 1.0, bottom_id)
bct_top  = DirichletBC(Q, 0.0, top_id)

# Write output files in VTK format:
u, p = z.split() # Do this first to extract individual velocity and pressure fields.
# Next rename for output:
u.rename("Velocity")
p.rename("Pressure")
# Create output file and select output_frequency:
output_file = File("output.pvd")
dump_period = 50
# Frequency of checkpoint files:
checkpoint_period = dump_period * 4
# Open file for logging diagnostic output:
f = open("params.log", "w")

# Now perform the time loop:
for timestep in range(0, max_timesteps):

    current_delta_t = delta_t
    if timestep != 0:
        delta_t.assign(compute_timestep(u, current_delta_t)) # Compute adaptive time-step
    time += float(delta_t)

    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0.
    solve(F_stokes==0, z, solver_parameters=stokes_solver_parameters, nullspace=Z_nullspace, transpose_nullspace=Z_nullspace, near_nullspace=Z_near_nullspace)

    # Temperature system:
    solve(F_energy==0, T_new, bcs=[bct_base,bct_top], solver_parameters=temperature_solver_parameters)

    # Write output:
    if timestep == 0 or timestep % dump_period == 0:
        output_file.write(u, p, T_new)

    # Compute diagnostics:
    u_rms               = sqrt(assemble(dot(u,u) * dx)) * sqrt(1./domain_volume)
    f_ratio             = rmin/rmax
    top_scaling         = -1.3290170684486309 # log(f_ratio) / (1.- f_ratio)
    bot_scaling         = -0.7303607313096079 # (f_ratio * log(f_ratio)) / (1.- f_ratio)
    nusselt_number_top  = (assemble(dot(grad(T_new),n) * ds_t) / assemble(Constant(1.0,domain=mesh)*ds_t) ) * top_scaling
    nusselt_number_base = (assemble(dot(grad(T_new),n) * ds_b) / assemble(Constant(1.0,domain=mesh)*ds_b) ) * bot_scaling
    energy_conservation = abs(abs(nusselt_number_top) - abs(nusselt_number_base))
    average_temperature = assemble(T_new * dx) / domain_volume

    # Calculate L2-norm of change in temperature:
    maxchange = sqrt(assemble((T_new - T_old)**2 * dx))

    # Log diagnostics:
    log_params(f, f"{timestep} {time} {maxchange} {u_rms} "
               f"{nusselt_number_base} {nusselt_number_top} "
               f"{energy_conservation} {average_temperature} ")

    # Leave if steady-state has been achieved:
    if maxchange < steady_state_tolerance:
        log("Steady-state achieved -- exiting time-step loop")
        break

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Checkpointing:
    if timestep % checkpoint_period == 0:
	# Checkpointing during simulation:
        checkpoint_data = DumbCheckpoint(f"Temperature_State_{timestep}", mode=FILE_CREATE)
        checkpoint_data.store(T_new, name="Temperature")
        checkpoint_data.close()

        checkpoint_data = DumbCheckpoint(f"Stokes_State_{timestep}", mode=FILE_CREATE)
        checkpoint_data.store(z, name="Stokes")
        checkpoint_data.close()

f.close()

# Write final state:
final_checkpoint_data = DumbCheckpoint("Final_Temperature_State", mode=FILE_CREATE)
final_checkpoint_data.store(T_new, name="Temperature")
final_checkpoint_data.close()

final_checkpoint_data = DumbCheckpoint("Final_Stokes_State", mode=FILE_CREATE)
final_checkpoint_data.store(z, name="Stokes")
final_checkpoint_data.close()

