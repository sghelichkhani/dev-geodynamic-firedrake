from firedrake import *
from mpi4py import MPI


x_max = 2.0; y_max = 1.0;
x_dis = 80; y_dis = 40;

# A built in periodic mesh 
mesh    = PeriodicRectangleMesh(x_dis, y_dis, x_max, y_max, direction='x')
X = x,y                = SpatialCoordinate(mesh)

# An extruded mesh with the same dimensions as the periodic box
mesh2d = IntervalMesh(ncells=x_dis, length_or_left=0, right=x_max) 
extruded_mesh = ExtrudedMesh(mesh2d, layers=y_dis, extrusion_type="uniform")
Qext   = FunctionSpace(extruded_mesh, "CG", 1)

# Function Space for periodic mesh 
Q       = FunctionSpace(mesh, "CG", 1) 

# The gaussian structure
blb_ctr_h = as_vector((1.5, 0.95))
blb_gaus = Constant(0.1)

# Set up temperature field and initialise based upon coordinates:
Tfield   = Function(Q, name="PeriodicTemperature")
Tfield.project((y_max - y) - 0.3*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2))

Tfieldext = Function(Qext, name="ExtrudedMesh")
Tfieldext.project(Tfield)

fi = File('Temperature.pvd')
fi.write(Tfieldext)

