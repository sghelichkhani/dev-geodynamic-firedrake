from firedrake import *

y_max = 1.0
x_max = 3.0

#mesh = PeriodicRectangleMesh(450, 150, x_max, y_max, direction='x')


# and Interval mesh of unit size 
mesh1d = IntervalMesh(100, length_or_left=0.0, right=x_max)

# extruding the base mesh "mesh1d" in the third dimension
mesh = ExtrudedMesh(mesh=mesh1d, layers=100, layer_height=y_max/100, extrusion_type='uniform', kernel=None, gdim=None)

Q2 = FunctionSpace(mesh, "CG", 2)

a_field = Function(Q2, name='field')

X = x, y = SpatialCoordinate(mesh)

a_field.interpolate(sqrt(x**2+y**2))

with CheckpointFile('a_checkpoint.h5', mode='w') as mesh_checkpoint:
    mesh_checkpoint.save_mesh(mesh)
    mesh_checkpoint.save_function(a_field, name='field')


# Now reload data
with CheckpointFile('a_checkpoint.h5', mode='r') as mesh_checkpoint:
    reloaded_mesh = mesh_checkpoint.load_mesh()

X_r = x_r, y_r = SpatialCoordinate(reloaded_mesh)

a_field_r = Function(Q2, name="field_reloaded")
a_field_r.interpolate(sqrt(x_r**2+y_r**2))

