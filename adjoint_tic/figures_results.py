import matplotlib
import numpy as np
import glob
import os
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm
import imageio

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
main_path = '/Users/sghelichkhani/Data/kerguelen/Workplace/dev-geodynamic-firedrake/adjoint_tic/rectangle_high_Ra/cases/'

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
     'BFGS_linesearch_1_SUPG_noplates': 'SUPG, top surface FS',\
     'BFGS_linesearch_1_SUPG_plates':   'SUPG, top surface imposed',\
     'BFGS_linesearch_1_noSUPG_plates': 'No SUP, top surface imposed',\
            }

ll_dirs = [os.path.join(main_path, i) for i in list(captions.keys())] 



iters    = [None]*len(ll_dirs) 
fvals    = [None]*len(ll_dirs) 
grads    = [None]*len(ll_dirs) 
nfvals   = [None]*len(ll_dirs) 
ngrads   = [None]*len(ll_dirs) 


# go through each and make sure there is an output file there
for i, case in enumerate(ll_dirs):
    if not os.path.isdir(case):
        continue
    try:
        trust_reg = 'Trust' in case
        iters[i], fvals[i], grads[i], nfvals[i], ngrads[i] =\
                    read_output_ROL(glob.glob(os.path.join(case,'*.out'))[0],\
                                   trust_reg_flg=trust_reg)
        fvals[i] /= np.max(fvals[i])
    except Exception:
        print('A problem with', case) 
        pass

def show_image(fi, ax):
    image =  imageio.imread(fi)
    ax.imshow(image)

y_0   = 0.39
x_0   = 0.00
del_x = 0.135*1.8
del_y = 0.3*1.8 
mar_y = 0.01
mar_x = 0.01
num_hor=6

for iteration in range(0, 51, 2):
    plt.close(1)
    fig = plt.figure(num=1, figsize=(24, 12))
    ax = fig.add_subplot(111, position=[x_0, y_0, del_x, del_y])
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    for i in ax.get_children():
        if isinstance(i, matplotlib.spines.Spine):
            i.set_color('red')
            i.set_linewidth(4)
    show_image('rectangle_high_Ra/PNGs/FWD_refmodel/fig_008.png', ax)
    bx = []
    for i, fi_name in enumerate(list(captions.keys())):
        axis_position = [x_0+(np.mod(i+1, num_hor-1))*(del_x+ mar_x), y_0-int((i+1)/(num_hor-1))*(del_y+mar_y), del_x, del_y]
        bx.append(fig.add_subplot(111, position=axis_position))
        try:
            show_image(os.path.join('rectangle_high_Ra/PNGs',fi_name,str('fig_fin_%2.2i.png') %ngrads[i][abs(np.array(iters[i])-iteration+1).argmin()]), bx[-1])
        except Exception:
            print('ouch')
            pass
        bx[-1].text(0.50, 0.95, captions[list(captions.keys())[i]], rotation=0, style='italic', ha='center', va='center', fontsize=18, transform=bx[-1].transAxes,\
                bbox={'facecolor':(0.9, 0.9, 0.9), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})
        bx[-1].set_axis_off()
    bx[-1].text(0.50, 0.95, str('Iteration %3.1i' %iteration), rotation=0, style='italic', ha='center', va='center', fontsize=32, transform=fig.transFigure,
                bbox={'facecolor':(0.9, 0.9, 0.9), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})
    fig_name=str('rectangle_high_Ra/PNGs/final_iter%3.3i.png'  %iteration)
    print(fig_name)
    fig.savefig(fig_name, dpi=100, bbox_inches='tight', pad_inches=0.1)
