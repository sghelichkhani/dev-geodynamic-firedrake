#!/usr/local/bin/pydoc3.9
import lib_optoutput
import matplotlib.pyplot as plt
import os.path
import glob
from matplotlib import cm
import numpy as np

scipy_fi_name = './'
scipy_misfit = lib_optoutput.readoutput(fi_name=scipy_fi_name, tag=2)

#
def on_pick(event):
    artist = event.artist
    xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
    x, y = artist.get_xdata(), artist.get_ydata()
    ind = event.ind
    print('Artist picked:', event.artist.get_label())
    print('x, y of mouse: {:.2f},{:.2f}'.format(xmouse, ymouse))
    print('Data point:', x[ind[0]], y[ind[0]])
    print(ind[0])
    print()

# What aire all the simunames
#simnames = [str('test_%3.3i' %i)for i in range(0, 96, 1)
simnames = glob.glob('test_*')
prms = [] 

for sim in simnames:
    with open(file=glob.glob(os.path.join(sim, sim+'.o*'))[0], mode='r') as fi:
        prms.append(fi.readline().split()[-4:])
prms = np.array(prms)

# set_conditions here
#mycondition = prms[:,0] == "Interpolation"
mycondition = [True]  * len(simnames)
temp_simnames = np.array(simnames)[mycondition]
temp_prms = prms[mycondition]


# go through them all
all_misfits = [lib_optoutput.readoutput(fi_name=glob.glob(os.path.join(sim, sim+'.o*'))[0], tag=0) for sim in temp_simnames]

plt.close(2)
fig = plt.figure(num=2, figsize=(25.6 , 12.91))
ax = fig.add_subplot(111)
ax.set_position([0.05, 0.05, 0.85, 0.85])
k= 0 
for i, misfit in enumerate(all_misfits):
    try:
        ax.plot(misfit[0][:,1], misfit[0][:,0], marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_misfits))(i), label=str("%s ," %temp_simnames[i])+" ,".join(np.array(temp_prms[i]))) 
        #ax.plot(range(len(misfit[0][:,1])), misfit[0][:,0], marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_misfits))(i), label=str("test_%3.3i" %i)+", ".join(np.array(temp_prms[i]))) 
    except:
        pass 

ax.plot(scipy_misfit[0][:,1], scipy_misfit[0][:,0], marker='o', markersize=4, picker=10, color='k') 
ax.set_ylim(1e-6,3e-3)
#ax.set_ylim(1e-6,3e-4)
ax.set_yscale('log')
#ax.legend(loc=(0.8, 0.0), ncol=2)
ax.grid()
fig.canvas.callbacks.connect('pick_event', on_pick)
fig.show()
