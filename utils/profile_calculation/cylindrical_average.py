from firedrake import *

rmin=1.22; rmax=2.22;

# Building the extruded mesh
mesh1d = CircleManifoldMesh(256*3, radius=rmin) 
mesh = ExtrudedMesh(mesh1d, layers=3*16, extrusion_type="radial")

# my coordinates
X = x,y                = SpatialCoordinate(mesh)

# Set up function space 
V       = FunctionSpace(mesh, "CG", 2) 

# The gaussian structure
blb_ctr_h = 1.6*as_vector((0.8, 0.5))
blb_gaus = Constant(0.1)

# Set up temperature field and initialise based upon coordinates:
T_field   = Function(V, name="OldTemperature")
T_field.project(0.5 + 0.5*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2))

# Sent by Steph:
#   mesh = V.mesh()
#   hcell, vcell = mesh.ufl_cell().sub_cells()
#   hele, _ = V.ufl_element().sub_elements()
#   vele = FiniteElement("R", vcell, 0)
#   ele = TensorProductElement(hele, vele)
#   V_1layer = FunctionSpace(mesh, ele)
# Modified to get the horizontal variation fixed
mesh = V.mesh()
hcell, vcell = mesh.ufl_cell().sub_cells()
hele, _ = V.ufl_element().sub_elements()
vele = FiniteElement("R", vcell, 0)
ele = TensorProductElement(hele, vele)
V_1layer = FunctionSpace(mesh, ele)

# func_ave is fixed horizontally
func_ave = Function(V_1layer)
func_ave.project(T_field)

# output
fi = File('Temperature.pvd')
fi.write(T_field, func_ave)

