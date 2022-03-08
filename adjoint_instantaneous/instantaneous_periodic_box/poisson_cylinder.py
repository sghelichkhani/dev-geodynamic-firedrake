"""
     This is a test for the poisson equation for gravity
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

# Geometric Constants:

r_max = 6370e3 
r_min = r_max - 2890e3

disc_order = 1
# setting of refinement level
refinement = 2**disc_order

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'

# Defining a single layer
m = CircleManifoldMesh(refinement*16*8, radius=r_min)

def radial_gaussion(x, mu, sigma):
    return 1/(sigma*numpy.sqrt(2*numpy.pi))*np.exp(-(x-mu)**2/(2*sigma**2))
def variable_radial_heights(nlayers):
    my_func = 1/(radial_gaussion(np.linspace(0, 1, nlayers), mu=0.0, sigma=0.22) +\
              radial_gaussion(np.linspace(0, 1, nlayers), mu=1.0, sigma=0.22))
    return np.ones(nlayers)*my_func/np.sum(my_func)*(r_max-r_min)

nlayers = refinement*16
rshl = np.ones(nlayers)/np.sum(np.ones(nlayers))*(r_max-r_min)

# Our extruded mesh
mesh = ExtrudedMesh(m, layers = rshl.size, layer_height=rshl, extrusion_type='radial')

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
    'ksp_type': 'gmres',
}


# Set up function spaces 
W    = FunctionSpace(mesh, "CG", 2)
Wvec    = VectorFunctionSpace(mesh, "CG", 1)

# Set up function space and associated test functions:
v       = TestFunction(W)   # test function for poisson equation
u       = TrialFunction(W)  # Trial Function for the poisson Equation
rho0    = Constant(4500)
G       = Constant(6.67430e-11) # Newtonian Gravitational Constant

## In case we want to have a spatially varying field
#rho = Function(W, name='density')
#rho.interpolate(interpolate(conditional(le(r, r_max), rho0, 0.0), W) * interpolate(conditional(ge(r, r_min-1000), rho0, 0.0), W))

# UFL: poisson equation
a = inner(grad(u), grad(v))*dx
L = 4 * G * pi * rho0 * v *dx

# Defining the analytical solution, page 187-188 of Thornton and Marion 
anl_sol = Function(W, name='AnalyticalSolution')
anl_sol.interpolate(-2*rho0*G*(r_max**2/2 - r_min**3/(3*r) - r**2/6))

# With us staying withing the boundaries we can use the analytical solution for the BSc 
bcs_poisson = DirichletBC(W, anl_sol, (top_id, bottom_id))
# If I decide the extend the boundaries, maybe I should change to these outside of the shell solutions
#bcu_poisson_inf = DirichletBC(W, -4/3*pi*rho0*G/r_max*(r_max**3-r_min**3), (top_id))
#bcu_poisson_0   = DirichletBC(W, -2*pi*rho0*G*(r_max**2-r_min**2), (bottom_id))

# solution
u = Function(W, name='potential')

# Solving for gravity
solve(a==L, u, bcs=[bcs_poisson], nullspace=VectorSpaceBasis(constant=True), solver_parameters=poisson_solver_parameters)

#
misfit = Function(W, name='Misfit')
misfit.assign(u - anl_sol)

# Write output of gravity 
g_file = File('out.pvd')
g_file.write(u, anl_sol, misfit)

