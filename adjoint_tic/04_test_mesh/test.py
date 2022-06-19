from firedrake import *

# dimensions of the mesh
y_max = 1.0
x_max = 3.0

# mesh generation
mesh = PeriodicRectangleMesh(450, 150, x_max, y_max, direction='x')

# function space 
Q2 = FunctionSpace(mesh, "CG", 2)
# function to map the distance from (0,0)
a_field = Function(Q2, name='field')

# UFL representation of coordinates
X = x, y = SpatialCoordinate(mesh)

# interpolating radius to function
a_field.interpolate(sqrt(x**2+y**2))

# writing out
with CheckpointFile('a_checkpoint.h5', mode='w') as mesh_checkpoint:
    mesh_checkpoint.save_mesh(mesh)
    mesh_checkpoint.save_function(a_field, name='field')


# Now reload data
with CheckpointFile('a_checkpoint.h5', mode='r') as mesh_checkpoint:
    reloaded_mesh = mesh_checkpoint.load_mesh()

# function space 
Q2_r = FunctionSpace(reloaded_mesh, "CG", 2)

X_r = x_r, y_r = SpatialCoordinate(reloaded_mesh)
#x_r, y_r = X_r[0], X_r[1]


a_field_r = Function(Q2_r, name="field_reloaded")
a_field_r.interpolate(sqrt(x_r**2+y_r**2))


pvd_file = File('test.pvd')
pvd_file.write(a_field)

pvd_file_r = File('test_r.pvd')
pvd_file_r.write(a_field_r)

