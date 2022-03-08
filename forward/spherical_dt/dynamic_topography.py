import gc 
from firedrake import *
from mpi4py import MPI
import math, numpy,  sys
from firedrake.petsc import PETSc

#########################################################################################################
################################## Some important constants etc...: #####################################
#########################################################################################################


# Dimensionalisation Constants:
mu_0    = 1e22                  # representative viscosity
L_0     = 2890e3                # rep. length
L_epsilon = 1e-2
u_0     = 1e-2/(3600*24*365)    # rep. velocity
p_0     = mu_0*u_0/L_0          # rep. pressure
Ra      = Constant(L_0**2/(mu_0*u_0))
delta_rho_um  = 4000                  # upper mantle density [kg/m^3]
delta_rho_cmb = 4000                  # upper mantle density [kg/m^3]

# Geometric Constants:
rmin = 1.22; rmax = 2.22;

# Projection solver parameters for nullspaces:
project_solver_parameters = {
    "snes_type": "ksponly",
    "ksp_type": "gmres",
    "pc_type": "sor",
    "mat_type": "aij",
    "ksp_rtol": 1e-12,
}

# Stokes Equation Solver Parameters:
stokes_solver_parameters = {
    'mat_type': 'matfree',
    'snes_type': 'ksponly',
    'ksp_type': 'preonly',
    'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_type': 'full',
    'fieldsplit_0': {
        'ksp_type': 'cg',
        'pc_type': 'python',
        'pc_python_type': 'firedrake.AssembledPC',
        'assembled_pc_type': 'gamg',
        'assembled_pc_gamg_threshold_scale': 1.0,
        'assembled_pc_gamg_threshold': 0.01,
        'assembled_pc_gamg_coarse_eq_limit': 800,
        # Note the following 5 options are required only when there is a velocity null space (which is the case here)
        "assembled_mg_coarse_ksp_type": "preonly",
        "assembled_mg_coarse_pc_type": "sor",
        "assembled_mg_coarse_ksp_max_it": 10,
        "assembled_mg_coarse_ksp_rtol": 0., # we always want 10 iterations
        "assembled_mg_coarse_ksp_atol": 0., # we always want 10 iterations
        'ksp_rtol': '1e-7',
        'ksp_test_transpose_null_space': None,
        'ksp_test_null_space': None,
        'mat_null_space_test_view': None,
        'ksp_monitor': None,        
        'ksp_converged_reason': None,
    },
    'fieldsplit_1': {
        'ksp_type': 'fgmres',
        'ksp_test_transpose_null_space': None,
        'ksp_test_null_space': None,
        'mat_null_space_test_view': None,
        'ksp_monitor': None,                'ksp_converged_reason': None,
        'pc_type': 'python',
        'pc_python_type': 'firedrake.MassInvPC',
        'Mp_ksp_type': 'cg',
        'Mp_pc_type': 'sor',
        'ksp_rtol': '1e-6',
    }
}


# Converting non-dimensional coordinates to dimensional
def dimensionalise(input_coords):
    return L_0*input_coords - 45.8e3


def RadialMu(rad):
    import numpy as np
    from scipy.interpolate import interp1d
    dpth, visc = np.loadtxt('./S10.abs', unpack=True)
    return interp1d(dpth*1e3, 10**visc/mu_0, fill_value="extrapolate")((rmax - rad)*L_0)

# Chosing boundary condition not on that boundary
class InteriorBC(DirichletBC):
    """DirichletBC applied to anywhere that is *not* on the specified boundary"""
    @utils.cached_property
    def nodes(self):
        return numpy.array(list(set(range(self._function_space.node_count)) - set(super().nodes)))


#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    PETSc.Sys.Print(*args)

def log_params(f, str):
    f.write(str + "\n")
    f.flush()

class rhs_data:
    def __init__(self, read_T=True, read_rho=True):
        import h5py
        from scipy.spatial import cKDTree
        self.read_rho = read_rho
        self.read_T   = read_T
        global u_xyz, data_x, data_y, data_z, dists, inds 
        data_fi     = h5py.File('./GS11_SL19_PY100_SO0_2900km.h5', 'r') 

        self.tree = cKDTree(np.column_stack((np.array(data_fi.get('coords/x')),\
                                             np.array(data_fi.get('coords/y')),\
                                             np.array(data_fi.get('coords/z')))))
        if read_T:  
            self.data_temperature    = np.array(data_fi.get('full_fields/temperature'))
        if read_rho:
            self.data_rho       = np.array(data_fi.get('full_fields/density'))
        data_fi.close()

    def __call__(self, mesh_coords, nns=10):
        dists, inds = self.tree.query(dimensionalise(mesh_coords), k=nns)
        dists[dists[:,0]<1e-9,1:] = 1e2; dists[dists[:,0]<1e-9,0] = 1.0
        return np.sum((1/dists)*self.data_temperature[inds], axis=1)/np.sum((1/dists),axis=1),\
                np.sum((1/dists)*self.data_rho[inds], axis=1)/np.sum((1/dists),axis=1)

        

#def model(ref_level,radial_layers):    
ref_level = 3
radial_layers = 32
# Mesh and associated physical boundary IDs:
mesh2d = IcosahedralSphereMesh(rmin, refinement_level=ref_level, degree=2)
mesh = ExtrudedMesh(mesh2d, radial_layers, (rmax-rmin)/radial_layers, extrusion_type='radial')
bottom_id, top_id = 'bottom', 'top'

# Stokes related constants:
X = x,y,z   = SpatialCoordinate(mesh)
r           = sqrt(X[0]**2+X[1]**2+X[2]**2)
theta       = atan_2(X[1], X[0]) # Theta (longitude - different symbol to Zhong)
phi         = atan_2(sqrt(X[0]**2+X[1]**2), X[2])  # Phi (co-latitude - different symbol to Zhong)
rhat        = as_vector((X[0]/r, X[1]/r, X[2]/r)) # Radial unit vector (in direction opposite to gravity)
n           = FacetNormal(mesh)

#########################################################################################################
################################## Geometry and Spatial Discretization: #################################
#########################################################################################################

# Set up function spaces - currently using the P2P1 element pair :
V       = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
Vscl    = FunctionSpace(mesh, "CG", 2)       # Velocity function space (vector)
W       = FunctionSpace(mesh, "CG", 1)       # Pressure function space (scalar)
Wvec    = VectorFunctionSpace(mesh, "CG", 1) # VectorFunction Space for spatial coordinates

# Output the number of vertices that we have
log('Number of vertices:',  W.dim())

# Set up mixed function space and associated test functions:
Z       = MixedFunctionSpace([V, W])
N, M    = TestFunctions(Z)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z       = Function(Z) # a field over the mixed function space Z.
u, p    = split(z) # can we nicely name mixed function space fields?

#########################################################################################################
############################################ Setup Equations ############################################
#########################################################################################################

# Equation in weak (ufl) form - note that continuity equation is added here - need to better understand why:
# Set up in residual form to ensure non-linear solvers are used.
g       = Constant(10.0)
C_ip    = Constant(20)
#mu      = Function(Vscl, name='Viscosity') # Constant viscosity
# coordinates in the same basis function as velocity
V_xyz   = interpolate(X, V)
V_r     = interpolate(r, Vscl)

# density
rhop = Function(Vscl, name='Density')

# rhs_data is a class that reads in the input data
rhs_model = rhs_data(read_T=True, read_rho=True)

# call for interpolation, giving in the coordinates, nns are the number of surrounding points used for interpolation
_, rhop.dat.data[:]  = rhs_model(V_xyz.dat.data[:], nns=5)
#mu.dat.data[:] =  RadialMu(V_r.dat.data[:])

mu = Constant(1.0)

# get rid of the class to free up some memory 
del rhs_model; gc.collect()

def tau(u):     return mu * (grad(u)+transpose(grad(u))) # Strain-rate tensor:
def trac(u,p):  return dot(tau(u),n) - p*n
def trac_field(u,p):  return dot(tau(u), conditional(gt(r, rmin+L_epsilon), -rhat, +rhat)) - p*conditional(gt(r, rmin+L_epsilon), -rhat, +rhat)

# nitsche free slip BCs
nitsche_fs  = - dot(N,n)*dot(n,trac(u,p))*ds_tb - dot(u,n)*dot(n,trac(N,M))*ds_tb\
        + C_ip*((2+1)*(2+3)/3)*FacetArea(mesh)/CellVolume(mesh)*dot(u,n)*dot(N,n)*ds_tb 

# Three possible rotation axes in 3 dimensions
# around x
x_rotZ = Function(Z)
x_rotV, _ = x_rotZ.split()
x_rotV.interpolate(as_vector((0, X[2], -X[1])))

# around y
y_rotZ = Function(Z)
y_rotV, _ = y_rotZ.split()
y_rotV.interpolate(as_vector((-X[2], 0, X[0])))

# around z
z_rotZ = Function(Z)
z_rotV, _ = z_rotZ.split()
z_rotV.interpolate(as_vector((-X[1], X[0], 0)))

# constant nullspace for pressure 
constW = VectorSpaceBasis(constant=True)

# Generating the nullspace basis and normalizing it
nullspaceV = VectorSpaceBasis([x_rotV, y_rotV, z_rotV])
nullspaceV.orthonormalize()

# Setting the mixed nullspace
nullspaceZ = MixedVectorSpaceBasis(Z, [nullspaceV, constW])

# Sorting out near nullspaces 
# three cartesian directions, x , y, z
# and three rotations
nns_rx = Function(V)
nns_ry = Function(V)
nns_rz = Function(V)
nns_x = Function(V)
nns_y = Function(V)
nns_z = Function(V)

nns_rx.interpolate(as_vector((    0, X[2], -X[1])))
nns_ry.interpolate(as_vector((-X[2],    0,  X[0])))
nns_rz.interpolate(as_vector((-X[1], X[0],     0)))
nns_x.interpolate(Constant([1., 0., 0.]))
nns_y.interpolate(Constant([0., 1., 0.]))
nns_z.interpolate(Constant([0., 0., 1.]))

# Constructing the basis
nns_VectorSpaceBasis = VectorSpaceBasis([nns_rx, nns_ry, nns_rz, nns_x, nns_y, nns_z])
nns_VectorSpaceBasis.orthonormalize()

# Construction of the nearnullspace for Z
nns_MixedVectorSpace = MixedVectorSpaceBasis(Z, [nns_VectorSpaceBasis, Z.sub(1)])

# Finalise equations
F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
F_stokes += -g * Ra * rhop * dot(N, rhat) * dx
F_stokes += - div(u)* M * dx
F_stokes += nitsche_fs

# Write output files in VTK format:
u, p = z.split() # Do this first to extract individual velocity and pressure fields:
u.rename('Velocity')
p.rename('Pressure')

u_file = File('velocity.pvd')
p_file = File('pressure.pvd')

# Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
# and the RHS == 0.
solve(F_stokes==0, z, solver_parameters=stokes_solver_parameters,  appctx={'mu': mu}, nullspace=nullspaceZ, transpose_nullspace=nullspaceZ, near_nullspace=nns_MixedVectorSpace)

## dimensionalising fields
u.project(u*u_0, solver_parameters=project_solver_parameters)
p.project(p*p_0, solver_parameters=project_solver_parameters)

# Computing dynamic topography
tau_rr = Function(Vscl, name='Dynamic Topography')
tau_rr.project(dot(trac_field(u, p), rhat)/(g*conditional(gt(r, rmin+L_epsilon), delta_rho_um, delta_rho_cmb)), solver_parameters=project_solver_parameters)

# domain  DirichletBC applied to anywhere that is *not* on the specified boundary 
tau_bc = InteriorBC(Vscl, 0.0, [top_id, bottom_id])
tau_bc.apply(tau_rr)

# for a p2 functionspace, we have nlayers*2+1 layers 
radial_tau = mesh.comm.allreduce(np.array([np.average(tau_rr.dat.data[i::radial_layers*2+1]) for i in range(radial_layers*2+1)]), op=MPI.SUM)/mesh.comm.size

# Map it onto our T_ave vector
tau_ave = Function(Vscl, name="average_tau")
tau_ave.dat.data[:] = np.array([[radial_tau[i] for i in range(radial_layers*2+1)]
                                for j in range(int(tau_ave.dat.data.shape[0]/(radial_layers*2+1)))]).reshape(np.shape(tau_ave.dat.data[:]))

# Remove average dynamic topography field 
tau_rr.project(tau_rr - tau_ave, solver_parameters=project_solver_parameters)

# Write output:
res_fi =  File('results.pvd')
res_fi.write(u, p, tau_rr,  rhop, tau_ave)

theta_field = Function(Vscl).project(theta, solver_parameters=project_solver_parameters)
phi_field   = Function(Vscl).project(phi, solver_parameters=project_solver_parameters)
r           = Function(Vscl).project(r, solver_parameters=project_solver_parameters)

checkpoint_data = DumbCheckpoint("Results", mode=FILE_CREATE)
checkpoint_data.store(theta_field, name="Azimuth")
checkpoint_data.store(phi_field, name="Longitude")
checkpoint_data.store(r, name="Radius")
checkpoint_data.store(tau_rr, name="DynamicTopography")
checkpoint_data.close()
