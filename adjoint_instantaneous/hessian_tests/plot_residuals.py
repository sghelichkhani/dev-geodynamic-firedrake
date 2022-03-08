

def read_optimization(fi_name):
    with open(fi_name, mode='r') as fi:
        k = 0
        fval = []
        iters = []
        nfcalc = []
        ndcalc = []
        real_final_state = []
        regular_term = []
        for line in fi:
            if str('%i' %k) in line[0:5]:
                try: 
                    fval.append(float(line.split()[1]))
                    iters.append(int(line[0:7]))
                    if k!=0:
                        nfcalc.append(int(line.split()[-3]))
                        ndcalc.append(int(line.split()[-2]))
                    else:
                        nfcalc.append(int(0))
                        ndcalc.append(int(0))
                    k = k +1 
                except:
                    pass
            if 'final state:' in line:
                real_final_state.append(float(line[len('final state:')+1:]))
            if 'regularisation:' in line:
                regular_term.append(float(line[len('regularisation:')+1:]))
    return iters, fval, nfcalc, ndcalc, real_final_state, regular_term

import glob 

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
plt.close(1)
fig = plt.figure(num=1, figsize=(14,7))
ax = fig.add_subplot(111, position=[0.10, 0.1, 0.6, 0.8])
all_simus = glob.glob('./hessian_17_MARCH/*.out')
all_simus.sort()
max_xaxis=50
for idx, i in enumerate(all_simus):
    try:
        beta = i.split('_')[-3]
    except:
        beta= 'No regularization'
    if i.split('_')[-2] == '0':
        linestyle='-'
        type_txt = 'BFGS'
    else:
        linestyle = '--'
        type_txt = 'BFGS with $H_0$'
    try:
        memory = i.split('_')[-1].split('.')[0]
    except:
        raise ValueError('something wrong with the naming') 
        
    if memory == '50':
        marker = 'o'
    else: 
        marker = '*'
    myiters, fvals, nfcalc, ndcalc, real_final_state, regular_term = read_optimization(fi_name=i)
    ax.plot(np.array(nfcalc)+np.array(ndcalc), fvals, label=str(r'%s, $\beta=%s$ memory: %s' %(type_txt, beta, memory)), linestyle=linestyle, color=cm.get_cmap('jet', len(all_simus))(idx), marker=marker)

#
ax.grid()
ax.set_yscale('log')
ax.legend(ncol=1, loc =[1.00, 0.0], fontsize=14, framealpha=1.0)
ax.set_xticks(list(np.arange(0, 500, 50)))
ax.set_xlim([-0.2, 300])
#ax.set_ylabel(r'$\sqrt{\frac{1}{2}(T(t=t_f)-T_{REF})^2\, dx}$')
ax.set_ylabel(r'$\sqrt{\frac{1}{2}(T(t=t_f)-T_{REF})^2\, dx} + \sqrt{\frac{1}{2}\nabla(t=t_I)\cdot\nabla(t=t_I)\, dx}$')
ax.set_xlabel('functional and derivative calc [\#]')
fig.show()
fig.savefig('./fevals_numfd.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
#fig.savefig('./real_misfit_evals.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
