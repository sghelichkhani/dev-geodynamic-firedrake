import lib_read_in
import matplotlib.pyplot as plt


fi_name = '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa_DefROL/cases/BFGS_linesearchCubic/BFGS_linesearchCubic.out'

misfit, fwd, der = lib_read_in.readoutput(fi_name, 1)
#
#plt.close(1)
#fig = plt.figure(num=1)
#ax = fig.add_subplot(111)
#ax.plot(misfit[:,1], misfit[:,0])
#ax.set_yscale('log')
#fig.show()
