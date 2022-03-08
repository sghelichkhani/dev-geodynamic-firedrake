"""
   Instantenous flow 
"""

from firedrake import *
from mpi4py import MPI
import math, numpy
from firedrake.petsc import PETSc

def main():
    #logging.set_log_level(1)
    #logging.set_level(1)
    
    # Geometric Constants:
    r_max = 6370e3 
    r_min = r_max - 2890e3

    disc_order = 1

    # setting of refinement level
    refinement = 2**disc_order

    # Top and bottom ids, for extruded mesh
    top_id, bottom_id = 'top', 'bottom'
    #top_id, bottom_id = 2, 1 


    # Defining a single layer
    m = CircleManifoldMesh(refinement*16*8, radius=r_min)

    # Our extruded mesh
    mesh = ExtrudedMesh(m, layers = refinement*16, layer_height=(r_max-r_min)/(refinement*16), extrusion_type='radial')


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
 

    # Stokes Equation Solver Parameters:
    stokes_solver_parameters = {
        'snes_type': 'ksponly',
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'pc_factor_mat_solver_type': 'mumps',
        'mat_type': 'aij'
    }
    
    # Set up function spaces - currently using the P2P1 element pair :
    V    = VectorFunctionSpace(mesh, "CG", 2) # Velocity function space (vector)
    W    = FunctionSpace(mesh, "CG", 1) # Pressure and Gravitational Potential function space (scalar)
    Wvec = VectorFunctionSpace(mesh, "CG", 1) # Vectorial function for gravity 
    
    # Set up mixed function space and associated test functions:
    Z       = MixedFunctionSpace([V, W])
    N, M    = TestFunctions(Z)
    Q       = TestFunction(W) # test function for poisson equation
    
    # Set up fields on these function spaces - split into each component so that they are easily accessible:
    z    = Function(Z)  # a field over the mixed function space Z.
    u, p = split(z)     # can we nicely name mixed function space fields?
    q    = Function(W)
    
    # Density and RHS definitions for the Stokes problem 
    rho0    = Constant(4500)
    gr0     = Constant(10.0)
    alpha   = Constant(2e-5)
    G       = Constant(6.67430e-11)

    # Having a single hot blob on 1.5, 0.0
    blb_ctr_h = as_vector((0., 0.5*(r_max+r_min))) 
    blb_gaus = Constant(0.1*(r_max-r_min))
    
    # set up the driving forces: right hand side 
    rhs    = Function(W, name="rhs")
    
    # A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
    rhs.interpolate(rho0 * gr0 *(1  - alpha * 500*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2)))

    # output denisty 
    rhs_file = File('rhs.pvd')
    rhs_file.write(rhs)

    # Setting boundary condition at the bottom
    bcu_rhs = DirichletBC(W, rho0 * gr0, (top_id, bottom_id))
    bcu_rhs.apply(rhs)

    # weak Stokes in ufl form
    mu        = Constant(1e22) # Constant viscosity
    beta      = Constant(20.0) # Magnitude of the penalty term, Used in weak imposition of the free-slip condition 

    # deviatoric stresses
    def tau(u): return  mu * (grad(u)+transpose(grad(u)))
    
    # traction field
    def trac(u,p): return dot(tau(u),n) - p*n

    F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
    F_stokes += (dot(N,rhat)*rhs) * dx 
    F_stokes += - div(u)* M * dx
    F_stokes += -dot(N,n)*dot(n,trac(u,p))*ds_tb -dot(u,n)*dot(n,trac(N,M))*ds_tb +beta/h*dot(u,n)*dot(N,n)*ds_tb

    # definition of the poisson equation for gravity
    F_poisson = inner(grad(q), grad(Q)) *dx - 4*pi*G*rho0*Q * dx - gr0 * Q * ds_b


    # For some reason this only works here!!!
    u, p    = z.split() 
    u.rename('Velocity') 
    p.rename('Pressure')
    
    # Printing out the degrees of freedom 
    if mesh.comm.rank == 0:
        log('Number of nodes:', W.dim())
    
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    solve(F_stokes==0, z, solver_parameters=stokes_solver_parameters)

    # Solving for gravity
    solve(F_poisson==0, q, solver_parameters=poisson_solver_parameters)
    
    # Taking the nullspace of velocity out
    rot = as_vector((-y,x))
    coeff =  assemble(dot(u, rot) * dx(domain=mesh)) / assemble(dot(rot, rot)*dx(domain=mesh))
    u.project(u - rot*coeff)

    # We know pressure at the surface is zero, we don't have negative pressure
    p.project(p-mesh.comm.allreduce(p.dat.data.min(), MPI.MIN))

    # Computing normal stresses acting on the surface
    tau_kk = Function(W, name='NormStress')
    tau_kk.interpolate(dot(dot(tau(u), rhat), -rhat)/(rho0*gr0))
    
    # domain condition for normal stresses
    tau_bc = InteriorBC(W, 0.0, [top_id, bottom_id])
    tau_bc.apply(tau_kk)

    # computing gravity from potential equation
    gravity = Function(Wvec, name='gravity')
    gravity.interpolate(-grad(q))

    # Write output of gravity 
    g_file = File('gravity.pvd')
    g_file.write(gravity)


    # Write output files in VTK format:
    u_file = File('velocity.pvd')
    p_file = File('pressure.pvd')
    tau_file = File('taukk.pvd')
    u_file.write(u)
    p_file.write(p)
    tau_file.write(tau_kk)
   

    ## Generating the reference temperature field for the adjoint
    #checkpoint_data = DumbCheckpoint("final_state", single_file=True, mode=FILE_CREATE, comm=mesh.comm)
    #checkpoint_data.store(T_new)
    #checkpoint_data.close()


#### Print function to ensure log output is only written on processor zero (if running in parallel) ####
def log(*args):
    PETSc.Sys.Print(*args) 

# Chosing boundary condition not on that boundary
class InteriorBC(DirichletBC):
    """DirichletBC applied to anywhere that is *not* on the specified boundary"""
    @utils.cached_property
    def nodes(self):
        return numpy.array(list(set(range(self._function_space.node_count)) - set(super().nodes)))


if __name__=='__main__':
    main()

