"""
    Adjoint Reconstruction - Using the classic way inputing parameters, instead of definiting methods for ROL.Algirithm() 
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc
from firedrake_adjoint import *
from pyadjoint import MinimizationProblem, ROLSolver
from pyadjoint.tape import no_annotations, Tape, set_working_tape
from pyadjoint.optimization.rol_solver import ROLObjective, ROLVector
from pyadjoint.optimization.optimization_solver import OptimizationSolver 
import ROL as ROL
import time; 
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
steady_state_tolerance = 1e-7
max_num_timesteps      = 50
target_cfl_no          = 2.5
max_timestep           = 1.00

# Stokes related constants:
Ra                     = Constant(1e6)   # Rayleigh Number

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

#########################################################################################################
############################ T advection diffusion equation Prerequisites: ##############################
#########################################################################################################

# Final state, which will be used as reference for minimization, loaded from a file 
final_state = Function(Q, name='RefTemperature')
final_state_file = DumbCheckpoint("../../final_state", mode=FILE_READ)
final_state_file.load(final_state, 'Temperature')
final_state_file.close()

# Initial condition
T_ic   = Function(Q, name="T_IC")
# Let's start with the final condition
T_ic.project(final_state)

# Set up temperature field and initialise based upon coordinates:
T_old    = Function(Q, name="OldTemperature")

# Having a single hot blob on 1.5, 0.0
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


### Temperature, advection-diffusion equation
F_energy = Y * ((T_new - T_old) / delta_t) * dx + Y*dot(u,grad(T_theta)) * dx + dot(grad(Y),kappa*grad(T_theta)) * dx

# For some reason this only works here!
u, p    = z.split() # Do this first to extract individual velocity, pressure and lagrange multplier fields:
u.rename('Velocity') 
p.rename('Pressure')

# A simulation time to track how far we are
simu_time = 0.0

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

    # updating time
    simu_time += float(delta_t)

    # Set T_old = T_new - assign the values of T_new to T_old
    T_old.assign(T_new)

    # Updating Temperature
    log("Timestep Number: ", timestep, " Timestep: ", float(delta_t))

## Initialise functional
functional = assemble(0.5*(T_new - final_state)**2 * dx)

# Below are callbacks allowing us to access various field information (accessed through reducedfunctional).
class OptimisationOutputCallbackPost:
    def __init__(self):
        self.iter_idx = 0
        self.opt_file             = File('visual/opt_file.pvd') 
        self.T_ic_true            = Function(Q, name="InitTemperature_Ref")
        self.T_ic_copy            = Function(Q, name="InitTemperature")
        self.T_tc_copy            = Function(Q, name="FinTemperature")

        # Having a single hot blob on 1.5, 0.0
        blb_ctr_h = as_vector((0.5, 0.85)) 
        blb_gaus = Constant(0.04)
        
        # A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
        self.T_ic_true.interpolate(0.5 - 0.3*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2))


    def __call__(self, cb_functional, dj, controls):
        # output current control (temperature initial condition)
        self.T_ic_copy.assign(controls)
        # output current final state temperature
        self.T_tc_copy.assign(T_new.block_variable.checkpoint)
        
        #  Write out the fields
        self.opt_file.write(self.T_ic_copy, self.T_tc_copy)
        func_val = assemble((self.T_tc_copy-final_state)**2 * dx) 
        init_func_val = assemble((self.T_ic_true-self.T_ic_copy)**2 * dx) 
        #reg_val  = assemble(inner(grad(self.T_ic_copy-T_mean), grad(self.T_ic_copy - T_mean)) * dx) 
        #grad_val = assemble((self.grad_copy)**2 * dx)

        log(f'# Der {self.iter_idx}, ||Misfit|| {func_val}, ||Misfit IC|| {init_func_val}')
        self.iter_idx += 1

class ForwardCallbackPost:
   def __init__(self):
      self.fwd_idx = 0 
   def __call__(self, func_value, controls):
      self.fwd_idx +=1 
      log(f'# fwd {self.fwd_idx}, ||Misfit|| {func_value}')

# Initiate classes for the callbacks
local_cb_post = OptimisationOutputCallbackPost()
eval_cb_post = ForwardCallbackPost()

# Defining the object for pyadjoint
reduced_functional = ReducedFunctional(functional, control, eval_cb_post=eval_cb_post, derivative_cb_post=local_cb_post)

# Set up bounds, which will later be used to enforce boundary conditions in inversion:
T_lb     = Function(Q, name="LB_Temperature")
T_ub     = Function(Q, name="UB_Temperature")
T_lb.assign(0.2)
T_ub.assign(0.5)

### Optimise using ROL - note when doing Taylor test this can be turned off:
minp = MinimizationProblem(reduced_functional, bounds=(T_lb, T_ub))

# This is the classic way
params = {
        'General': {
              'Print Verbosity': 1,
              'Secant': {'Type': 'Limited-Memory BFGS', 'Maximum Storage': 10},
                    },
        'Step': {
           'Type': 'Line Search',
           'Line Search': {
                'Descent Method': {'Type': 'Quasi-Newton Method'},
                'Line-Search Method': {
                                'Type':  "Brent's", #'Cubic Interpolation ''Backtracking',#'Bisection',
                                'Backtracking Rate': 0.5,
                                'Bracketing Tolerance': 1.e-1,
                                'Bisection': {
                                    'Tolerance': 1e-1,
                                    'Iteration Limit': 20,
                                            },
                                        },
                                "Brent's": {
                                    'Tolerance': 1e0,
                                    'Iteration Limit': 3,
                                    'Run Test Upon Initialization': False,
                                            },
                'Curvature Condition': {
                                'Type': 'Strong Wolfe Conditions',
                                'General Parameter': 0.9,
                                'Generalized Wolfe Parameter': 0.6,
                                        },
                'Function Evaluation Limit': 20,
                'Sufficient Decrease Tolerance': 1e1,
                'Use Previous Step Length as Initial Guess': False,
                            },
                },
        'Status Test': {
            'Gradient Tolerance': 0,
            'Iteration Limit': 10,
                        }
        }

# overwritting ROLObjective to have a cache
class myROLObjective(ROLObjective):
    def __init__(self, rf, scale=1.0, f_cachesize=4, g_cachesize=2):
        super().__init__(rf, scale=scale)
        
        # cache size for functionals and gradients
        self.f_cachesize = f_cachesize
        self.g_cachesize = g_cachesize

        # cache for x, given for functional and gradient calculations 
        self.fx = []
        self.gx = []

        # cache for result of functional and gradient calculations 
        self.fvals       = []
        self.grads = []

        # Sia: to see the actual gradient that is being ued (not l2)
        self.gradFile    = File(filename='./gradients/gradient.pvd')
        self.g_pvd_field = Function(Q, name="gradient")

    # functional value is accessed here 
    def value(self, x, tol):
        return self.val 

    # updating the gradient g.dat
    def gradient(self, g, x, tol):

        # check if x is already stored in cache
        idx = self.cachescan(self.gx, x)

        if idx==None:
            # in case the last forward run used a different x, rerun again
            if self.cachescan(self.fx, x) != 0:
                self.rf(x.dat)

            # cache x.dat 
            self.gx.insert(0, [Function(f.function_space()).assign(f) for f in x.dat])

            # gradient calculation
            init_time = time.perf_counter()
            super().gradient(g, x, tol)
            log(f"Elapsed time for grad calc {time.perf_counter() - init_time} sec")

            # cache g.dat 
            self.grads.insert(0, [Function(g.function_space()).assign(g) for g in g.dat])
            
            ## Write out recently computed gradient
            self.gradFile.write(self.g_pvd_field.assign(g.dat[0]))

        # if x is found in cache
        else:
            # idx is the index of gradient field in cache 
            [g.dat[i].assign(cg) for i, cg in enumerate(self.grads[idx])]

        # size control for cache
        while len(self.gx) >self.g_cachesize:
            self.gx.pop()
            self.grads.pop()

    # updating self.val which is passed to ROL by self.value 
    def update(self, x, flag, iteration):

        # check if x is already stored in cache 
        idx = self.cachescan(self.fx, x)

        if idx==None: 
            # cache x.dat 
            self.fx.insert(0, [Function(f.function_space()).assign(f) for f in x.dat])

            # update self.val
            init_time = time.perf_counter()
            super().update(x, flag, iteration)
            log(f"Elapsed time for func eval {time.perf_counter() - init_time} sec")
            
            # Storing fval 
            self.fvals.insert(0, self.val)
        else:
            # update control to the cache value 
            for i, value in enumerate(self.fx[idx]):
                self.rf.controls[i].update(value)
            # Update value 
            self.val = self.fvals[idx]

        # size control for cache
        while len(self.gx) > self.f_cachesize:
            self.fx.pop()
            self.fvals.pop()


    # scanning cache in case we already have the fields 
    def cachescan(self, cache, iterx):

        # if cache is empty 
        if len(cache)==0:
            return None 

        # idx is the index of the field in cache 
        idx = None

        # check if we already have x
        for j, xc in enumerate(cache):
            if np.sum([float(assemble((xc[i] - iterx.dat[i])**2*dx)) for i, _ in enumerate(iterx.dat)]) <= numpy.finfo(float).eps:
                idx = j
                break
        return idx 


class myROLSolver(ROLSolver):
    def __init__(self, problem, parameters, inner_product="L2"):
       """
       Generate a ROL solver that uses myROLObective instead of ROLObjective

       The argument inner_product specifies the inner product to be used for
       the control space.

       """

       OptimizationSolver.__init__(self, problem, parameters)
       self.rolobjective = myROLObjective(problem.reduced_functional)
       x = [p.tape_value() for p in self.problem.reduced_functional.controls]
       self.rolvector = ROLVector(x, inner_product=inner_product)
       self.params_dict = parameters

       self.bounds = self._ROLSolver__get_bounds()
       self.constraints = self._ROLSolver__get_constraints()


with stop_annotating():
    # set up ROL problem
    rol_solver = myROLSolver(minp, params, inner_product="L2")
    sol = rol_solver.solve()

    # Save the optimal temperature_ic field 
    ckpt_T_ic = DumbCheckpoint("T_ic_optimal",\
            single_file=True, mode=FILE_CREATE,\
                               comm=mesh.comm)
    ckpt_T_ic.store(sol)
    ckpt_T_ic.close()
