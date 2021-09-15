import numpy as np
import glob
import os
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm

def read_output_ROL(fi_name, trust_reg_flg=False):
    """
        Read the result of an adjoint optimisation with Rapid Optimisation Library
    """
    i = 0
    with open(file=fi_name, mode='r') as fi:
        # A list of all the items
        myiters = []; myfvals = []; mygrad  = []; mynfval  = []; myngrad  = []; 
    
        for line in fi:
            line_spl = line.split()
            try:
                if '{}'.format(i) == line_spl[0]:
                    try:
                        if trust_reg_flg:
                            nfval_idx=5
                            max_zeroth=4
                        else:
                            nfval_idx=4
                            max_zeroth=3
                        mynfval.append(float(line_spl[nfval_idx]))
                        myngrad.append(float(line_spl[nfval_idx+1]))
                    except Exception:
                        if len(line_spl)==max_zeroth:
                            mynfval.append(0)
                            myngrad.append(0)
                    myiters.append(i)
                    myfvals.append(float(line_spl[1]))
                    mygrad.append(float(line_spl[2]))
                    i+=1
            except Exception:
                pass 
    
    return myiters, myfvals, mygrad, mynfval, myngrad


# where all the different cases are located
main_path = '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_high_Ra/cases'

# Captions
captions = {
     #'w_reg_i_BFGS':            'BFGS Cubical Default prms',\
     #'w_reg_02_BFGS':           'BFGS Cubical $10^{-3}$ Regularisation',\
     #'w_reg_03_BFGS':           'BFGS Cubical $10^{-4}$ Regularisation',\
     #'NLCG_Fletcher_Rieves':    'NLCG Fletcher Rieves',\
     #'NLCG_Polak_Ribiere':      'NLCG Polak Ribiere',\
     #'SteepestDescent_Backtracking_StrongWolf' : 'Steepest Descent Backtracking Strong Wolf',\
     #'SteepestDescent_Backtracking_Wolf': 'Steepest Descent Backtracking Wolf',\
     #'BFGS_linesearh_1':        'BFGS Cubical MaxFeval=5',\
     #'BFGS_H0_LineSearch1':     'BFGS Cubical 20 Max Eval $H_0$ ',\
     #'BFGS_linesearh_2':        'BFGS Cubical MaxFeval=2',\
     #'BFGS_linesearh_3':        'BFGS Cubical tol=$10^{-4}$',\
     #'BFGS_linesearh_4':        'BFGS Cubical tol=$10^{-3}$',\
     #'BFGS_linesearh_5':        'BFGS Cubical use last step',\
     #'BFGS_H0_LineSearch5':     'BFGS Cubical $H_0$ use last step',\
     #'BFGS_Trust_Region':       'BFGS Trust Region',\
     #'BFGS_H0_Trust_Region':    'BFGS Trust Region $H_0$',\
      'BFGS_linesearch_1_SUPG_noplates': 'BFGS with Cubic interpolation, SUPG, top surface FS',\
      'BFGS_linesearch_1_SUPG_plates': 'BFGS with Cubic interpolation, SUPG, top surface imposed',\
      'BFGS_linesearch_1_noSUPG_plates': 'BFGS with Cubic interpolation, No SUP, top surface imposed',\
            }


## get the name of all the cases
#ll_dirs = glob.glob(os.path.join(main_path, '*'))
#ll_dirs.sort()

ll_dirs = [os.path.join(main_path, i) for i in list(captions.keys())] 



plt.close(1)
fig = plt.figure(num=1, figsize=(18, 8))
ax = fig.add_subplot(111, position=[0.1, 0.1, 0.6, 0.7])

# go through each and make sure there is an output file there
for i, case in enumerate(ll_dirs):
    if not os.path.isdir(case):
        continue
    try:
        trust_reg = 'Trust' in case
        iters, fvals, grads, nfvals, ngrads =\
                    read_output_ROL(glob.glob(os.path.join(case,'*.out'))[0],\
                                   trust_reg_flg=trust_reg)
        fvals = fvals/np.max(fvals)
        #ax.plot(np.array(nfvals)+np.array(ngrads), fvals, label=case[case.rfind('/')+1:].replace('_', ' '),\
        #                    color=cm.get_cmap('jet', len(ll_dirs))(i))
        #if np.max(iters[-1])>50:
        #    ax.plot(np.array(nfvals[:50])+np.array(ngrads[:50]), fvals[:50], label=captions[list(captions.keys())[i]],\
        #                    color=cm.get_cmap('jet', len(ll_dirs))(i))
        #else:
        ax.plot(np.array(nfvals)+np.array(ngrads), fvals, label=captions[list(captions.keys())[i]],\
                        color=cm.get_cmap('jet', len(ll_dirs))(i))

    except Exception:
        print('A problem with', case) 
        pass
ax.set_xlim(0, 100)
ax.set_ylabel(r'Normalized $f$')
ax.text(0.50, 1.08, 'Up to 500 Iterations', rotation=0, style='italic', ha='center', va='center', fontsize=32, transform=ax.transAxes,\
        bbox={'facecolor':(0.9, 0.9, 0.9), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})
ax.set_xlabel(r'\# of $f$+ $\nabla g$ evals')
ax.legend(loc=(0.20, 0.7),framealpha=1.0, facecolor='w')
ax.grid()
ax.set_yscale('log')
#fig.show()
fig.savefig('./high_ra_box.png', dpi=100, bbox_inches='tight', pad_inches=0.1) 
#
#
