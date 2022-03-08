from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

m = UnitSquareMesh(5, 5)
mesh = ExtrudedMesh(m, 5, layer_height=[0.1, 0.1, 0.4, 0.1, 0.1], extrusion_type='uniform')
X = SpatialCoordinate(mesh)
V = VectorFunctionSpace(mesh, 'CG', 1)
out = interpolate(X, V)
fi = File('out.pvd')
fi.write(out)

