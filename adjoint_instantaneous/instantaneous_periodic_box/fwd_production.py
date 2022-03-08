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
    y_max = 2890e3
    x_max = 2890e3
    
    #  how many intervals along x/y directions 
    disc_n = 100
    
    # Top and bottom ids, for extruded mesh
    #top_id, bottom_id = 'top', 'bottom'
    top_id, bottom_id = 2, 1 
    
    # The mesh
    mesh = utility_meshes.PeriodicRectangleMesh(nx=disc_n, ny=disc_n,\
                                                Lx=x_max, Ly=y_max, direction='x',\
                                                quadrilateral=True, reorder=None,\
                                                distribution_parameters=None, diagonal=None)

    # setting up spatial coordinates
    X  =  x, y = SpatialCoordinate(mesh)
    
    # a measure of cell size
    h	  = sqrt(CellVolume(mesh))
    
    # setting up vertical direction
    yhat  = as_vector((0,1))
    
    # Stokes related constants:
    Ra                     = Constant(1e5)   # Rayleigh Number
    

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
    Q       = TestFunction(W)
    
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
    blb_ctr_h = as_vector((0.5*y_max, 0.5*y_max)) 
    blb_gaus = Constant(0.05*y_max)
    
    # set up the driving forces: right hand side 
    rhs    = Function(W, name="rhs")
    
    # A linear temperature profile from the surface to the CMB, with a gaussian blob somewhere
    rhs.interpolate(rho0 * gr0 *(1  - alpha * 500*exp(-0.5*((X-blb_ctr_h)/blb_gaus)**2)))
    
    # Setting boundary condition at the bottom
    bcu_rhs = DirichletBC(W, rho0 * gr0, (top_id, bottom_id))
    bcu_rhs.apply(rhs)
    # weak Stokes in ufl form
    mu        = Constant(1e22) # Constant viscosity
    
    # deviatoric stresses
    def tau(u): return  mu * (grad(u)+transpose(grad(u)))
    
    F_stokes  = inner(grad(N), tau(u)) * dx - div(N)*p * dx 
    F_stokes += (dot(N,yhat)*rhs) * dx 
    F_stokes += - div(u)* M * dx

    # definition of the poisson equation for gravity
    F_poisson = + inner(grad(q), grad(Q)) *dx - 4*pi*G*rho0*Q * dx - gr0 * Q * ds(bottom_id)

    # Setting free-slip BC for top and bottom
    bcu_topbase= DirichletBC(Z.sub(0).sub(1), 0.0, (top_id, bottom_id))
    
    # For some reason this only works here!!!
    u, p    = z.split() 
    u.rename('Velocity') 
    p.rename('Pressure')
    
    # Printing out the degrees of freedom 
    if mesh.comm.rank == 0:
        log('Number of nodes:', W.dim())
    
    # Solve system - configured for solving non-linear systems, where everything is on the LHS (as above)
    # and the RHS == 0. 
    solve(F_stokes==0, z, bcs=[bcu_topbase], solver_parameters=stokes_solver_parameters)

    # Solving for gravity
    solve(F_poisson==0, q, solver_parameters=poisson_solver_parameters)
    
    # Taking the nullspace of velocity out
    e_x = as_vector((1,0))
    coeff =  assemble(dot(u, e_x) * dx(domain=mesh)) / assemble(dot(e_x, e_x)*dx(domain=mesh))
    u.project(u - e_x*coeff)

    # We know pressure at the surface is zero, we don't have negative pressure
    p.project(p-mesh.comm.allreduce(p.dat.data.min(), MPI.MIN))

    # Computing normal stresses acting on the surface
    tau_kk = Function(W, name='NormStress')
    tau_kk.interpolate(dot(dot(tau(u), yhat), -yhat)/(rho0*gr0))
    
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
    t_file = File('rhs.pvd')
    tau_file = File('taukk.pvd')
    u_file.write(u)
    p_file.write(p)
    t_file.write(rhs)
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

