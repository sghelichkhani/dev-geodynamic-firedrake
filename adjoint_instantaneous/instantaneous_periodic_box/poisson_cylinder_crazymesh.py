"""
   Instantenous flow 
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

# Geometric Constants:

rcenter = 1.0
r_max = 6370e3 
r_min = r_max - 2890e3
r_inf = r_max * 10
disc_order = 1
# setting of refinement level
refinement = 2**disc_order

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'
#top_id, bottom_id = 2, 1 

# Defining a single layer
m = CircleManifoldMesh(refinement*16*8, radius=rcenter)

def radial_gaussion(x, mu, sigma):
    return 1/(sigma*numpy.sqrt(2*numpy.pi))*np.exp(-(x-mu)**2/(2*sigma**2))
def variable_radial_heights(nlayers):
    my_func = 1/(radial_gaussion(np.linspace(0, 1, nlayers), mu=0.0, sigma=0.22) +\
              radial_gaussion(np.linspace(0, 1, nlayers), mu=1.0, sigma=0.22))
    return np.ones(nlayers)*my_func/np.sum(my_func)*(r_max-r_min)
#rshl =  variable_radial_heights(refinement*16)

rshl = np.hstack((np.ones(4) * ((r_min - rcenter)/4), variable_radial_heights(refinement*16), np.ones(4)*((r_inf - r_max)/4)))

# Our extruded mesh
mesh = ExtrudedMesh(m, layers = rshl.size,\
                    layer_height=rshl,\
                    extrusion_type='radial')


# setting up spatial coordinates
X  =  x, y = SpatialCoordinate(mesh)

# a measure of cell size
h	  = sqrt(CellVolume(mesh))
n     = FacetNormal(mesh) 
# setting up vertical direction
r     = sqrt(x**2 + y**2)
rhat  = as_vector((x, y)) / r




# Poisson Equation Solve parameters 
poisson_solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_type': 'mumps',
    'mat_type': 'aij'
}

#solver_parameters={'ksp_type': 'cg'}


# Set up function spaces 
W    = FunctionSpace(mesh, "DG", 2)
#Wvec    = VectorFunctionSpace(mesh, "DG", 2)

# Set up function space and associated test functions:
v       = TestFunction(W) # test function for poisson equation
u       = TrialFunction(W)
rho0    = Constant(4500)
G       = Constant(6.67430e-11) # Newtonian Gravitational Constant
rho = conditional(le(r, r_max), rho0, 0.0)

# defining the poisson equation
poisson = inner(grad(u), grad(v))*dx  
bcu_poisson_inf = DirichletBC(W, 0, (top_id))
bcu_poisson_0   = DirichletBC(W, 100, (bottom_id))

# solution
u       = Function(W, name='potential')
# Solving for gravity
#solve(poisson==0, u, bcs=[bcu_poisson_0], solver_parameters=poisson_solver_parameters)

# Write output of gravity 
g_file = File('solutions.pvd')
g_file.write(u, interpolate(rho, W))

