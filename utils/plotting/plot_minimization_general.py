import matplotlib.pyplot as plt
import os

# Definig a class for all the output files
class optimisation_output:
    def __init__(self, fi_name, tag, caption='  '):
        from lib_read_in import readoutput
        self.fi_name = fi_name
        self.tag = tag
        self.caption = caption
        self.misfit, self.fwd_calcs, self.adj_calcs = readoutput(self.fi_name, self.tag)

repos_path = os.path.join(os.getenv('HOME'), 'Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_extrude_examples')
fi_extension = 'cases/BFGS_linesearchCubic/BFGS_linesearchCubic.out'
sims = [\
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_single_proc_Scipy', fi_extension)        , tag = 2, caption= 'Scipy-l2'),
       ##optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_Scipy', fi_extension)             , tag = 2, caption= 'NonlinearSolve Scipy'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa', fi_extension)                   , tag = 0, caption= 'ROL-L2-Cubic Int'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_ROLl2', fi_extension)             , tag = 0, caption= 'ROL-l2-Cubic Int'),
       ##optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_ROL_Scaled', fi_extension)               , tag = 0, caption= 'ROL scaled with 4e3'),
       ##optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_ROL_Scaled_Higher', fi_extension) , tag = 0, caption='ROL scaled with 4e3, and l2'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_DefROL', fi_extension)            , tag = 1, caption= 'ROL Defaults (l2 Trust-Region)'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_Q2_HighRa_ROLl2', fi_extension)          , tag = 0, caption= 'T:Q2, ROL-l2-Cubic Int'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_Q2_HighRa_ROLL2', fi_extension)          , tag = 0, caption= 'T:Q2, ROL-L2-Cubic Int'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_DBD_HighRa_ROLL2', fi_extension)         , tag = 0, caption= 'T:Dirichlet BCs, ROL-L2-Cubic Int'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_RolManualScaled', fi_extension)   , tag = 0, caption= 'ROL-L2-Cubic Int, (Manually scaling gradient by 1e-4)'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'rectangle_extruded_ns_HighRa_RolManualScaled_Bisection', fi_extension), tag = 0, caption= 'ROL-L2-Bisection, (Manually scaling gradient by 1e-4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/Scipy_InitConsistent/Scipy_InitConsistent.out'),            tag = 2, caption= 'Scipy-L2-Backtracking, (Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/ROL_CubicInterpolation/ROL_CubicInterpolation.out'),        tag = 0, caption= 'ROL-L2-CubicInterpolation, (Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/ROL_Backtracking/ROL_Backtracking.out'),                    tag = 0, caption= 'ROL-L2-Backtracking (Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/ROL_Backtracking_FreeGuess/ROL_Backtracking_FreeGuess.out'),tag = 0, caption= 'ROL-L2-Backtracking with Free guess(Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/ROL_CubicInterpolation_lastalpha/ROL_CubicInterpolation_lastalpha.out'),tag = 0, caption= 'ROL-L2-Cubic Interpolation Last Alpha (Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test_bigger/cases/Scipy_BFGS/Scipy_BFGS.out'),                         tag = 2, caption= 'Bigger Blob - Scipy-L2 Backtracking (Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/Scipy_NoBounds/Scipy_NoBounds.out'),                        tag = 2, caption= 'Scipy-L2 Backtracking No Bounds(Scaling F by 1e+4)'),
       optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/Scipy_NoBounds_Saturation/Scipy_NoBounds_Saturation.out'),  tag = 2, caption= 'Scipy-L2 Backtracking No Bounds But Saturated'),
       #optimisation_output(fi_name=os.path.join(repos_path, 'test/cases/Manual/Manual.out'),                                    tag = 3, caption= 'Scipy-L2, (Scaling F by 1e+4)'),
       ]

plt.close(1)
fig = plt.figure(num=1, figsize=(19.16, 12.0))
x_0 = 0.1; y_0 = 0.65; del_x = 0.2; del_y = 0.2; mar_x= 0.05; mar_y=0.08; num_cols = 3
ax = []
for idx, output in enumerate(sims):
    ax.append(fig.add_subplot(111))
    ax[-1].set_position([x_0+(idx%num_cols)*(mar_x+del_x), y_0 - int(idx/num_cols)*(del_y+mar_y), del_x, del_y])
    ax[-1].plot(output.misfit[output.adj_calcs[:-3,2].astype('int'),1], output.adj_calcs[:-3,0], color='blue', label='Terminal Cond. (Functional)')
    ax[-1].plot(output.misfit[output.adj_calcs[:-3,2].astype('int'),1], output.adj_calcs[:-3,1], color='red', label='Initial Cond.')
    #ax[-1].plot(vv[:,0], vv[:,2], color='red', label='Initial Misfit')
    ax[-1].set_xlabel('Functional + Derivative [\#]')
    ax[-1].set_ylabel(r'$\frac{1}{2}\ (T - T_{ref})\ dx$')
    ax[-1].set_yscale('log')
    ax[-1].set_xlim([0, 1800])
    ax[-1].set_ylim([1e-7, 1e-1])
    ax[-1].grid()
    ax[-1].text(0.50, 1.05, str('%s' %output.caption), rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=ax[-1].transAxes,\
                bbox={'facecolor':(1.0, 1.0, 1.0), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})

ax[1].legend(loc=[0.5, 1.2])
#fig.show()
fig.savefig('./minimisation_results.png', dpi=100, bbox_inches='tight', pad_inches=0.1)
