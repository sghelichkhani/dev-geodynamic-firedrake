# This is not a complete code!
##################### For computing the average in an Extruded Mesh ##############################
ref_level              = 4
nlayers                = 32
# Mesh and associated physical boundary IDs:
mesh2d = IcosahedralSphereMesh(rmin, refinement_level=ref_level, degree=2)
mesh = ExtrudedMesh(mesh2d, nlayers, (rmax-rmin)/nlayers, extrusion_type='radial')

# for a p2 functionspace, we have nlayers*2+1 layers 
radial_temp = np.array([mesh.comm.allreduce(np.sum(rad.dat.data[i::nlayers*2+1])/(Q.dim()/(nlayers*2+1)), op=MPI.SUM) for i in range(nlayers*2+1)])

# Map it onto our T_ave vector
T_ave = Function(Q)
T_ave.dat.data[:] = np.array([[radial_temp[i] for i in range(nlayers*2+1)] for j in range(int(len(T_ave.dat.data[:])/(nlayers*2+1)))]).reshape(np.shape(T_ave.dat.data[:]))
#################################################################################################


