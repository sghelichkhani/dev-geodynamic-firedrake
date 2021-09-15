#!/usr/local/bin/pydoc3.9
import lib_optoutput
import matplotlib.pyplot as plt
import os.path
import glob
from matplotlib import cm
import numpy as np


class optresults():
    def __init__(self, fi_name, tag, caption):
        self.fi_name = fi_name
        self.tag     = tag
        self.caption = caption 
    def misfit(self):
        return lib_optoutput.readoutput(fi_name=self.fi_name, tag=self.tag)  
        
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
all_results= [
    #optresults(fi_name='./test_scipy/Scipy.o27393728', tag=2, caption="Scipy"),\
    optresults(fi_name='./test_scipy_scaledF/ScipyScaledF.o27393729', tag=2, caption="Scipy (Functional*1e4)"),\
    optresults(fi_name='./test_scipy_scaledg/ScipyScaledg.o27393730', tag=2, caption="Scipy (Gradient 1e-4)"),\
    optresults(fi_name='./test_TR_001/TrustRegion.o27394195', tag=1, caption="ROL BFGS \& TR "),\
    optresults(fi_name='./test_TR_002/TrustRegionFScaled.o27400356', tag=1, caption="ROL BFGS \& TR (Functional*1e4)"),\
    optresults(fi_name='./test_TR_003/test_TR_003.o27400402', tag=1, caption="ROL BFGS \& TR (Gradient*1e4)"),\
    ]

plt.close(2)
fig = plt.figure(num=2, figsize=(25.6 , 12.91))
ax = fig.add_subplot(111)
ax.set_position([0.05, 0.05, 0.85, 0.85])
k= 0 
for i, res in enumerate(all_results):
    to_plot = res.misfit()
    #ax.plot(to_plot[0][:200,1], to_plot[0][:,0]/to_plot[0][0,0],  marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_results))(i), label=str(res.caption)) 
    ax.plot(range(len(to_plot[0][:,0])), to_plot[0][:,0]/to_plot[0][0,0],  marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_results))(i), label=str(res.caption)) 
ax.set_ylim(1e-6,1.1)
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
fig.canvas.callbacks.connect('pick_event', on_pick)
fig.show()
