#!/usr/local/bin/python3.9
import lib_optoutput
import matplotlib.pyplot as plt
import os.path
import glob
from matplotlib import cm
import numpy as np
import itertools
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


class optresults(object):
    def __init__(self, fi_name):
        self.fi_name = fi_name
        self.tag     = self.get_tag() 
        # Get our caption
        self.get_caption()
    def misfit(self):
        return lib_optoutput.readoutput(fi_name=self.fi_name, tag=self.tag) 
    def get_tag(self):
        if "LS" in self.fi_name:
            return 0
        elif "TR" in self.fi_name or "SR" in self.fi_name:
            return 1
        elif "scipy" in self.fi_name:
            return 2
        else:
            return 1
            #raise ValueError(f"TAG is not classified for {self.fi_name}")
    def get_caption(self):
        with open(file=self.fi_name, mode='r')  as fi_in:
            temp = fi_in.readline()
            if (len(temp.split()) > 1 and "LS" in self.fi_name) or ("LINESEARCHMETHOD" in temp):
                self.caption = " ".join(temp.split(" ")[4:]).replace('\n','')
            else:
                self.caption = temp.replace('\n', '')
            #if "_scaledG_" in self.fi_name:
            #    self.caption += ",Gradient scaled with 1e-4"

            #if "_TR_" in self.fi_name:
            #    self.caption ="ROL TR with BFGS window:"+temp.split()[-1]
        self.caption = self.caption.replace('_','\_')




all_results  = [optresults(fi) for fi in glob.glob('*/*.out')]
plt.close(2)
fig = plt.figure(num=2, figsize=(8.0,  7.5))
ax = fig.add_subplot(111)
ax.set_position([0.20, 0.19, 0.60, 0.50])
k= 0 
for i, res in enumerate(all_results):
    to_plot = res.misfit()
    try:
        ax.plot(to_plot[0][:,1], to_plot[0][:,0]/np.max(to_plot[0][:,0]),  marker='o', markersize=3, picker=10,color='red', label='Final Misfit') 
        ax.plot(to_plot[0][:,1][to_plot[2][:,2].astype('int')], to_plot[2][:,1]/np.max(to_plot[2][:,1]),  marker='o', markersize=3, picker=10,color='blue', label='Initial Misfit') 
        #ax.plot(range(len(to_plot[0][:,1])), to_plot[0][:,0],  marker='o', markersize=3, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_results))(i), label=res.caption) 
    except:
        print(res.fi_name)
        pass
ax.set_ylabel('Normalised Misfit f')
ax.set_xlim([0,500])
ax.set_xlabel('Cost [Functional + Gradient Calculations]') 
ax.set_yscale('log')
ax.legend(loc=(0.15, 1.00), ncol=2)
ax.grid()
fig.canvas.callbacks.connect('pick_event', on_pick)
fig.show()
#fig.savefig('./comparisons_cost.png', dpi=200, bbox_inches='tight', pad_inches=0.1)
