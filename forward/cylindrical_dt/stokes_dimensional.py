from firedrake import *
from mpi4py import MPI
import math, numpy, sys
from firedrake.petsc import PETSc

# We define the principal dimensions, by which we non-dimensionalise our fields
L_nd = 2890e3   # depth of the mantle [km]
mu_nd = 1e20     # the minimum viscosity [Pa s]
U_nd = 1e-2/(3600*24*365) # representative velocity in the mantle: 1[cm/year]-> [m/s]
rho_um = 3200
# Set up geometry and key parameters:
rmin, rmax = 1.22, 2.22
epsilon = 1e-1
k          = Constant(2)  # radial degree
nn         = Constant(4)  # wave number (n is already used for FacetNormal)
g          = Constant(9.8)

# Stokes related constants (note that since these are included in UFL, they are wrapped inside Constant):
#mu        = Constant(1.0) # Constant viscosity
Ra        = Constant(L_nd**2/(mu_nd*U_nd))
C_ip      = Constant(20) # The fudge factor for interior penalty term used in weak imposition of BCs
dim       = 2

#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    PETSc.Sys.Print(*args)

#### File for logging output diagnostics through simulation.
def log_params(f, str):
    f.write(str + "\n")
    f.flush()


# Chosing boundary condition not on that boundary
class InteriorBC(DirichletBC):
    """DirichletBC applied to anywhere that is *not* on the specified boundary"""
    @utils.cached_property
    def nodes(self):
        return numpy.array(list(set(range(self._function_space.node_count)) - set(super().nodes)))

#########################################################################################################
######################################### Solver parameters:  ###########################################
#########################################################################################################

# Stokes Equation Solver Parameters:
stokes_solver_parameters = {
    "mat_type": "matfree",
    "snes_type": "ksponly",
    "ksp_type": "preonly",
    "pc_type": "fieldsplit",
    "pc_fieldsplit_type": "schur",
    "pc_fieldsplit_schur_type": "full",
    "fieldsplit_0": {
        "ksp_type": "cg",
        "pc_type": "python",
        "pc_python_type": "firedrake.AssembledPC",
        "assembled_pc_type": "gamg",
        "assembled_pc_gamg_threshold_scale": 1.0,
        "assembled_pc_gamg_threshold": 0.01,
        "assembled_pc_gamg_coarse_eq_limit": 800,
        # Note the following 5 options are required only when there is a velocity null space (which is the case here)
        "assembled_mg_coarse_ksp_type": "preonly",
        "assembled_mg_coarse_pc_type": "sor",
        "assembled_mg_coarse_ksp_max_it": 10,
        "assembled_mg_coarse_ksp_rtol": 0., # we always want 10 iterations
        "assembled_mg_coarse_ksp_atol": 0., # we always want 10 iterations
        "ksp_rtol": 1e-5,     
        "ksp_converged_reason": None,
    },
    "fieldsplit_1": {
        "ksp_type": "fgmres",
        "ksp_converged_reason": None,
        "pc_type": "python",
        "pc_python_type": "firedrake.MassInvPC",
        "Mp_ksp_type": "cg",
        "Mp_pc_type": "sor",
        "ksp_rtol": 1e-3,
    }
}

# Projection solver parameters for nullspaces:
project_solver_parameters = {
    "snes_type": "ksponly",
    "ksp_type": "gmres",
    "pc_type": "sor",
    "mat_type": "aij",
    "ksp_rtol": 1e-12,
}

disc_n = 4

# Construct a circle mesh and then extrude into a cylinder:
mesh1d = CircleManifoldMesh(disc_n*256, radius=rmin)
mesh   = ExtrudedMesh(mesh1d, layers=disc_n*16, extrusion_type="radial")
# Physical boundary IDs from extruded mesh:
bottom_id, top_id = "bottom", "top"

# Set-up P2 version (popped) of mesh
P1    = FunctionSpace(mesh, "CG", 1)
V     = VectorFunctionSpace(mesh, "CG", 2)
x, y  = SpatialCoordinate(mesh)
r     = sqrt(x**2 + y**2)
r_p1  = interpolate(r, P1)
xy_p2 = interpolate(as_vector((x,y))/r*r_p1, V)
mesh  = Mesh(xy_p2)

# now redefine everything based on the super/iso parametric P2 mesh
X = x, y  = SpatialCoordinate(mesh)
n     = FacetNormal(mesh)
h     = sqrt(CellVolume(mesh))
r     = sqrt(x**2 + y**2)
rhat  = as_vector((x, y)) / r
domain_volume = assemble(1.*dx(domain=mesh))

# For RHS:
phi       = atan_2(y, x)
rhop   = (4500 - 50*r**k/rmax**k*cos(nn*phi))

#########################################################################################################
################################## Geometry and Spatial Discretization: #################################
#########################################################################################################

# Set up function spaces - currently using the P2P1 element pair :
V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
W    = FunctionSpace(mesh, "CG", 1) # Pressure function space (scalar)
Wvec = VectorFunctionSpace(mesh, "CG", 1) # used for coordinates of pressure nodes

# Set up mixed function space and associated test functions:
Z    = MixedFunctionSpace([V, W])
N, M = TestFunctions(Z)

# Set up fields on these function spaces - split into each component so that they are easily accessible:
z      = Function(Z) # a field over the mixed function space Z.
u, p   = split(z)
u_, p_ = z.split()

#########################################################################################################
############################################ Setup Equations ############################################
#########################################################################################################

# Define a variable viscosity
mu = Constant(1.0) 

# deviatoric stresses
def tau(u): return  mu * (grad(u)+transpose(grad(u)))

# traction field 
def trac(u,p): return dot(tau(u),n) - p*n

def trac_new(u,p): return dot(tau(u),rhat) - p*rhat

# nitsche free slip BCs
# (p+1)*(p+d)/d; with p:= polynomial degree, d:= dimension => (2+1)(2+2)/2
nitsche_fs  = - dot(N,n)*dot(n,trac(u,p))*ds_tb - dot(u,n)*dot(n,trac(N,M))*ds_tb + C_ip*((2+1)*(2+dim)/dim)*FacetArea(mesh)/CellVolume(mesh)*dot(u,n)*dot(N,n)*ds_tb

# UFL for Stokes equations (note that given we are imposing free-slip BCs weakly, we have done some integration by parts:
F_stokes  = inner(grad(N), tau(u)) * dx -div(N) * p * dx
F_stokes += - g * Ra * rhop * dot(N, rhat) * dx
F_stokes += -div(u) * M * dx # Continuity equation
F_stokes += nitsche_fs

# Nullspaces and near-nullspaces:
rotZ     = Function(Z)
rotV, _  = rotZ.split()
rotV.interpolate(as_vector((-y, x)))
rotVSpace = VectorSpaceBasis([rotV])
rotVSpace.orthonormalize()

# Constant nullspace for pressure
p_nullspace = VectorSpaceBasis(constant=True)

# Setting the mixed nullspace
Z_nullspace = MixedVectorSpaceBasis(Z, [rotVSpace, p_nullspace])

# Generating near_nullspaces for GAMG, starting with constant modes:
c0V = Function(V)
c0V.interpolate(Constant([1., 0.]))
c1V = Function(V)
c1V.interpolate(Constant([0., 1.]))
# Rotation around Z axis (using rotV from above here messes with nullspaces due to orthonormalisation)
cr0V = Function(V)
cr0V.interpolate(as_vector((-y, x)))

V_near_nullspace = VectorSpaceBasis([cr0V, c0V, c1V])
V_near_nullspace.orthonormalize()
Z_near_nullspace = MixedVectorSpaceBasis(Z, [V_near_nullspace, Z.sub(1)])

# Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
# and the RHS == 0.
solve(F_stokes==0, z, solver_parameters=stokes_solver_parameters, appctx={"mu": mu}, nullspace=Z_nullspace, transpose_nullspace=Z_nullspace, near_nullspace=Z_near_nullspace)


# Computing dynamic topography
tau_rr = Function(W, name='NormalStress')
tau_rr.project(mu_nd*U_nd/L_nd/g/rho_um*dot(trac_new(u, p), conditional(gt(r, rmin+epsilon), -rhat, +rhat)), solver_parameters=project_solver_parameters)

# Applying the Interioir boundary condition
tau_bc = InteriorBC(W, 0.0, [top_id, bottom_id])
tau_bc.apply(tau_rr)


# for a p2 functionspace, we have nlayers*2+1 layers 
radial_tau = np.array([mesh.comm.allreduce(np.sum(tau_rr.dat.data[i::disc_n*16+1])/(W.dim()/(disc_n*16+1)), op=MPI.SUM) for i in range(disc_n*16+1)])
# Map it onto our T_ave vector
tau_ave = Function(W, name="average_tau")
tau_ave.dat.data[:] = np.array([[radial_tau[i] for i in range(disc_n*16+1)] for j in range(int(len(tau_ave.dat.data[:])/(disc_n*16+1)))]).reshape(np.shape(tau_ave.dat.data[:]))

tau_rr.project(tau_rr - tau_ave, solver_parameters=project_solver_parameters)


# projecting density on W
rhop_ = Function(W, name='Temperature')
rhop_.project(rhop, solver_parameters={'ksp_type':'gmres'})

# Write output files in VTK format:
u_.rename("Velocity")
p_.rename("Pressure")
out_file    = File("varmu_results.pvd")

u_dm = Function(V, name='Velocity')
u_dm.project(u_*U_nd)
p_dm = Function(W, name='Pressure')
p_dm.project(p_*mu_nd*U_nd/L_nd)

# Write output:
out_file.write(u_dm, p_dm, tau_rr, tau_ave, rhop_)

