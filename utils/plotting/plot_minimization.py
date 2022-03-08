import matplotlib.pyplot as plt
from lib_read_in import read_out_file 
import os

sims = [\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_single_proc_Scipy/',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa_Scipy',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa_ROLl2',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_ROL_Scaled',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa_ROL_Scaled_Higher',\
       '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples/rectangle_extruded_ns_HighRa_DefROL',\
       ]

caps = [\
       'NS all, Ra = 1e8, LinearSolve Scipy',\
       'NS all, Ra = 1e8, NonlinearSolve Scipy',\
       'NS all, Ra = 1e8, ROL L2 Cubic Interpolatation',\
       'NS all, Ra = 1e8, ROL l2 Cubic Interpolatation',\
       'NS all, Ra = 1e8, ROL scaled with 4e3',\
       'NS all, Ra = 1e8, ROL scaled with 4e3, and l2',\
       'NS all, Ra = 1e8, ROL Def Parameters',\
       ]
#fi_name'/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_Q2/cases/BFGS_linesearchCubic_wPlates/BFGS_linesearchCubic_wPlates.out'         = '/Users/sghelichkhani/Workplace/develope_G-ADOPT/adjoint_tic/rectangle_extruded/cases/BFGS_linesearchCubic_wPlates/out.dat'
#fi_name = '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_fs/cases/BFGS_linesearchCubic_wPlates/BFGS_linesearchCubic_wPlates.out'

plt.close(1)
fig = plt.figure(num=1, figsize=(19.16, 12.0))
x_0 = 0.1; y_0 = 0.65; del_x = 0.2; del_y = 0.2; mar_x= 0.05; mar_y=0.08; num_cols = 3
ax = []
for idx, fi_name in enumerate(sims):
    vv = read_out_file(    os.path.join(os.path.join(fi_name, 'cases/BFGS_linesearchCubic/BFGS_linesearchCubic.out')))
    ax.append(fig.add_subplot(111))
    ax[-1].set_position([x_0+(idx%num_cols)*(mar_x+del_x), y_0 - int(idx/num_cols)*(del_y+mar_y), del_x, del_y])
    ax[-1].plot(vv[:,0], vv[:,1], color='blue', label='Final Misfit (Functional)')
    ax[-1].plot(vv[:,0], vv[:,2], color='red', label='Initial Misfit')
    ax[-1].set_xlabel('Derivative Calls [\#]')
    ax[-1].set_ylabel('1/2 (T - T\_ref)')
    ax[-1].set_yscale('log')
    ax[-1].set_xlim([0, 1000])
    ax[-1].set_ylim([1e-9, 1e-2])
    ax[-1].grid()
    ax[-1].text(0.50, 1.05, str('%s' %caps[idx]), rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=ax[-1].transAxes,\
                bbox={'facecolor':(1.0, 1.0, 1.0), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})

ax[-1].legend()
fig.show()
#fig.savefig('./scipy_tools_Ra1e8_l2.png', dpi=100, bbox_inches='tight', pad_inches=0.1)
