"""
     This is a test for the poisson equation for gravity
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

# Geometric Constants:
r_max = 1.0 
r_min = 0.5 

# setting of refinement level
ref_level = 3  
nlayers = 2 ** ref_level 

# Top and bottom ids, for extruded mesh
top_id, bottom_id = 'top', 'bottom'

# Defining a single layer
m = IcosahedralSphereMesh(r_min, refinement_level=ref_level, degree=2)
# Our extruded mesh
mesh = ExtrudedMesh(m, layers = nlayers, layer_height=(r_max-r_min)/nlayers,\
                    extrusion_type='radial')

# setting up spatial coordinates
X  =  x, y, z = SpatialCoordinate(mesh)
r     = sqrt(x**2 + y**2 + z**2)

# Poisson Equation Solve parameters 
poisson_solver_parameters = {'snes_type': 'ksponly','ksp_type': 'gmres'}

# Set up function spaces 
W    = FunctionSpace(mesh, "CG", 2)

# Set up function space and associated test functions:
v       = TestFunction(W)   # test function for poisson equation
u       = TrialFunction(W)  # trial Function for the poisson Equation
rho0    = Constant(4500)    # density of the mantle
G       = Constant(6.67430e-11) # Newtonian Gravitational Constant

# UFL: poisson equation
a = inner(grad(u), grad(v))*dx
L = 4 * G * pi * rho0 * v *dx

# Defining the analytical solution, page 187-188 of Thornton and Marion 
anl_sol = Function(W, name='AnalyticalSolution')
anl_sol.interpolate(-4*pi*rho0*G*(r_max**2/2 - r_min**3/(3*r) - r**2/6))

# Staying withing the boundaries we use the analytical solution as Dirichlet BSc 
bcs_poisson = DirichletBC(W, anl_sol, (top_id, bottom_id))

# gravitational potential 
u = Function(W, name='potential')

# solving for gravity
solve(a==L, u, bcs=[bcs_poisson], solver_parameters=poisson_solver_parameters)

# computing our  misfit
misfit = Function(W, name='Misfit')
misfit.assign(u - anl_sol)

# Write output of gravity 
g_file = File('out.pvd')
g_file.write(u, anl_sol, misfit)

