from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake_adjoint import *
from firedrake.petsc import PETSc
from pyadjoint.optimization.optimization import minimize
from pyadjoint.tape import no_annotations, Tape, set_working_tape
import ROL
import sys

#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################

beta = Constant(float(sys.argv[1]))
init_hessian = int(sys.argv[2])

# Mesh and associated physical boundary IDs:
mesh    = PeriodicRectangleMesh(160, 80, 2, 1, direction='x')
bottom_id, top_id = 1, 2  # with PeriodicRectangleMesh: 1,2 for bottom and top

# Global Constants:
steady_state_tolerance = 1e-8
max_timesteps          = 2
model_time             = 0.0

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number
k                      = Constant((0,1))
X                      = SpatialCoordinate(mesh)

# Temperature related constants:
delta_t                = Constant(5.0e-6) # Initial time-step
kappa                  = Constant(1.0)  # Thermal diffusivity

# Temporal discretisation - Using a Crank-Nicholson scheme where theta_ts = 0.5:
theta_ts               = 0.5

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    if mesh.comm.rank == 0:
        PETSc.Sys.Print(*args)
log(init_hessian)

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

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Set up T_ic which will try to invert for in the adjoint simulation
T_ic     = Function(Q, name="Temperature_IC")
T_bg     = conditional(X[1] >= 0.96,  0.5 - ( (X[1]-0.96) / 0.04 )*0.5, 0.5)
T_bg_bc  = conditional(X[1] <= 0.04,  0.5 + ( (0.04-X[1]) / 0.04 )*0.5, T_bg)
T_final  = conditional(T_bg_bc < 0.0, 0.0, T_bg_bc)
T_ic.interpolate(T_final)

# Set up bounds, which will later be used to enforce boundary conditions in inversion:
T_lb     = Function(Q, name="LB_Temperature")
T_ub     = Function(Q, name="UB_Temperature")
T_lb.assign(0.0)
T_ub.assign(1.0)

# Define control for inversion:
control = Control(T_ic)

# Set up temperature field and initialise based upon coordinates:
T_old   = Function(Q, name="OldTemperature")
T_new   = Function(Q, name="Temperature")

# Temporal discretisation - Using a theta scheme:
T_theta = theta_ts * T_new + (1-theta_ts) * T_old

#########################################################################################################
################################ Load Checkpoint Info ###################################################
#########################################################################################################

final_state     = Function(Q, name="FinalState")
checkpoint_data = DumbCheckpoint("./Final_State",mode=FILE_READ)
checkpoint_data.load(final_state,name="Temperature")
checkpoint_data.close()

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
bcv_val_top  = Function(V).interpolate(as_vector((2000.0, 0.0)))
bcv_top      = DirichletBC(Z.sub(0), bcv_val_top, top_id)
#bcv_top      = DirichletBC(Z.sub(0).sub(1), 0.0, top_id)
bcv_val_base = Function(V).interpolate(as_vector((0.0, 0.0)))
bcv_base     = DirichletBC(Z.sub(0), bcv_val_base, bottom_id)

### Next deal with Temperature advection-diffusion equation: ###
delta_x      = sqrt(CellVolume(mesh))
Y_SUPG       = Y + dot(u,grad(Y)) * (delta_x / (2*sqrt(dot(u,u))))
F_energy     = Y_SUPG * ((T_new - T_old) / delta_t) * dx + Y_SUPG*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx
bct_base     = DirichletBC(Q, 1.0, bottom_id)
bct_top      = DirichletBC(Q, 0.0, top_id)

### Next set up boundary conditions for T_IC control. Without these, boundary values are not satisfied in the inversion:
bct_lb_base = DirichletBC(Q, 1.0, bottom_id)
bct_ub_top  = DirichletBC(Q, 0.0, top_id)

# Overwrite BC of control with correct values:
bct_lb_base.apply(T_lb)
bct_ub_top.apply(T_ub)

T_old.assign(T_ic)
T_new.assign(T_ic)

# Write output files in VTK format:
u, p = z.split() # Do this first to extract individual velocity and pressure fields:
u._name = 'Velocity' # This doesn't seem to work!
p._name = 'Pressure'
u_file     = File('velocity.pvd')
p_file     = File('pressure.pvd')
t_file     = File('temperature.pvd')
mu_file    = File('viscosity.pvd')

###################################################################################################################
#######################################    Adjoint Callbacks   ####################################################
###################################################################################################################

# Initialise functional
functional = 0.0

# Below are callbacks relating to the adjoint solutions (accessed through solve).
tape = get_working_tape()

##############################################################################################################

# Now perform the time loop:
for timestep in range(0,max_timesteps):

    log("Timestep Number: ", timestep, "Time: ", model_time)

    if(timestep % 10 == 0):
        # Write output:
        u_file.write(u)
        p_file.write(p)
        t_file.write(T_new)
        mu_field.interpolate(mu)
        mu_file.write(mu_field)

    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0.
    solve(F_stokes==0, z, bcs=[bcv_base,bcv_top], solver_parameters=stokes_solver_parameters)

    # Temperature system:
    solve(F_energy==0, T_new, bcs=[bct_base,bct_top], solver_parameters=temperature_solver_parameters)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Update model_time
    model_time = model_time + float(delta_t)

###################################################################################################################
########################################## Evaluate Functional ####################################################
###################################################################################################################

final_state_functional = assemble((T_new - final_state)**2 * dx)
log("Final State Functional:", final_state_functional)

log('Magnitude of regulariser:', float(beta))
regularisation_term = assemble(beta*dot(grad(T_ic),grad(T_ic)) * dx)
log("Regularisation Functional:", regularisation_term)

functional += final_state_functional
functional += regularisation_term
log("Overall Functional:", functional)

### Reduced funcational - functional with respect to the control.
reduced_functional = ReducedFunctional(functional, control)

### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))

class InitHessian(ROL.InitBFGS):
    def __init__(self, M, direct=True):
        ROL.InitBFGS.__init__(self, M)

        self.solver_params = {
            "ksp_type": "gmres",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "mat_type": "aij",
        }
        # K matrix
        self.M = dot(grad(Y), grad(TrialFunction(Q))) * dx

        # regular linear solver
        self.solver = LinearSolver(assemble(self.M), solver_parameters=self.solver_params)

    @no_annotations
    def applyH0(self, Hv, v):
        self.solver.solve(Hv.dat[0], v.dat[0])
        l = Hv.norm() / v.norm()
        Hv.scale(1 / l)
        self.scaleH0(Hv)

    @no_annotations
    def applyB0(self, Bv, v):
        Bv.dat[0] = assemble(dot(grad(Y), grad(v.dat[0])) * dx)
        l = Bv.norm() / v.norm()
        Bv.scale(1 / l)
        self.scaleB0(Bv)

class StatusReport(ROL.StatusTest):
    def __init__(self, params, vector):
        ROL.StatusTest.__init__(self, params)
        self.vector = vector
        self.T_copy              = Function(Q, name="Temperature")
        self.opt_t_final_file    = File('opt_temperature_final.pvd')
        self.opt_t_init_file     = File('opt_temperature_initial.pvd')

    @no_annotations
    def check(self, status):
        log("final state: {}".format(assemble(0.5*(T_new.block_variable.checkpoint - final_state)**2 * dx)))
        log("regularisation: {}".format(assemble(0.5*beta*dot(grad(self.vector.dat[0]), grad(self.vector.dat[0])) * dx)))

        # Writing out final condition
        self.T_copy.assign(T_new.block_variable.checkpoint)
        self.opt_t_final_file.write(self.T_copy)

        ## Writing out initial condition
        self.T_copy.assign(self.vector.dat[0])
        self.opt_t_init_file.write(self.T_copy) 
        return ROL.StatusTest.check(self, status)

# Orginial parameters by Steph
params = {
        'General': {
            'Secant': {'Type': 'Limited-Memory BFGS', 'Maximum Storage': 3},
            'Print Verbosity': 1},
        'Step': {
            'Line Search': {
                'User Defined Initial Step Size': True,
                'Step Size': 0.5 
                           }
                },
        'Status Test': {
            'Gradient Tolerance': 1e-12,
            'Iteration Limit': 20,
        }
    }

# set up ROL problem
rol_solver = ROLSolver(minp, params)
params = ROL.ParameterList(params, "Parameters")

if init_hessian == 1:
    log('Using user-defined initial hessian')
    rol_secant = InitHessian(3) # maximum storage
else:
    log('Using unit matrix for initial hessian')
    rol_secant = ROL.InitBFGS(3) # maximum storage
rol_step = ROL.TrustRegionStep(rol_secant, params)
#rol_step = ROL.LineSearchStep(params, rol_secant)
#rol_status = ROL.StatusTest(params)
rol_status = StatusReport(params, rol_solver.rolvector)
rol_algorithm = ROL.Algorithm(rol_step, rol_status)
log("functional before going for minimisation: {}".format(assemble(0.5*(T_new - final_state)**2 * dx)))

with stop_annotating():
    rol_algorithm.run(
    	rol_solver.rolvector,
        rol_solver.rolobjective,
        rol_solver.bounds # only if we have a bounded problem
            )

optimal_ic = rol_solver.problem.reduced_functional.controls.delist(rol_solver.rolvector.dat)

# Save the optimal temperature_ic field 
ckpt_T_ic = DumbCheckpoint("T_ic_optimal",\
        single_file=True, mode=FILE_CREATE,\
                           comm=mesh.comm)
ckpt_T_ic.store(optimal_ic)
ckpt_T_ic.close()


