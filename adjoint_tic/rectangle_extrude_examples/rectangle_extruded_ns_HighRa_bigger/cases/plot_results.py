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
        # Get our caption
        self.get_caption()
    def misfit(self):
        return lib_optoutput.readoutput(fi_name=self.fi_name, tag=self.tag)  
    def get_caption(self):
        with open(file=self.fi_name, mode='r')  as fi_in:
            temp = fi_in.readline()
            if len(temp.split()) > 1 and "LS" in self.fi_name:
                self.caption = " ".join(temp.split(" ")[4:]).replace('\n','')
            else:
                self.caption = temp.replace('\n', '')
        
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
    #optresults(fi_name="test_scipy_scaledF/test_scipy_scaledF.o27503231", tag=2, caption="40"),\
    optresults(fi_name="test_scipy_scaledg/test_scipy_scaledg.o27503180", tag=2, caption="41"),\
    optresults(fi_name="test_LS_001/test_LS_001.o27502540"              , tag=0, caption="0 "),\
    optresults(fi_name="test_LS_002/test_LS_002.o27502541"              , tag=0, caption="1 "),\
    optresults(fi_name="test_LS_003/test_LS_003.o27502542"              , tag=0, caption="2 "),\
    optresults(fi_name="test_LS_004/test_LS_004.o27502543"              , tag=0, caption="3 "),\
    optresults(fi_name="test_LS_005/test_LS_005.o27502544"              , tag=0, caption="4 "),\
    optresults(fi_name="test_LS_006/test_LS_006.o27502545"              , tag=0, caption="5 "),\
    optresults(fi_name="test_LS_007/test_LS_007.o27502546"              , tag=0, caption="6 "),\
    optresults(fi_name="test_LS_008/test_LS_008.o27502547"              , tag=0, caption="7 "),\
    optresults(fi_name="test_LS_009/test_LS_009.o27502548"              , tag=0, caption="8 "),\
    optresults(fi_name="test_LS_010/test_LS_010.o27502549"              , tag=0, caption="9 "),\
    optresults(fi_name="test_LS_011/test_LS_011.o27502550"              , tag=0, caption="10"),\
    #optresults(fi_name="test_LS_F1e4_000/test_LS_F1e4_000.o27502957"    , tag=0, caption="11"),\
    #optresults(fi_name="test_LS_F1e4_001/test_LS_F1e4_001.o27502958"    , tag=0, caption="12"),\
    #optresults(fi_name="test_LS_F1e4_002/test_LS_F1e4_002.o27502959"    , tag=0, caption="13"),\
    #optresults(fi_name="test_LS_F1e4_003/test_LS_F1e4_003.o27502960"    , tag=0, caption="14"),\
    #optresults(fi_name="test_LS_F1e4_004/test_LS_F1e4_004.o27502961"    , tag=0, caption="15"),\
    #optresults(fi_name="test_LS_F1e4_005/test_LS_F1e4_005.o27502962"    , tag=0, caption="16"),\
    #optresults(fi_name="test_LS_F1e4_006/test_LS_F1e4_006.o27502963"    , tag=0, caption="17"),\
    #optresults(fi_name="test_LS_F1e4_007/test_LS_F1e4_007.o27502964"    , tag=0, caption="18"),\
    #optresults(fi_name="test_LS_F1e4_008/test_LS_F1e4_008.o27502965"    , tag=0, caption="19"),\
    #optresults(fi_name="test_LS_F1e4_009/test_LS_F1e4_009.o27502966"    , tag=0, caption="20"),\
    #optresults(fi_name="test_LS_F1e4_010/test_LS_F1e4_010.o27502967"    , tag=0, caption="21"),\
    #optresults(fi_name="test_LS_F1e4_011/test_LS_F1e4_011.o27502968"    , tag=0, caption="22"),\
    optresults(fi_name="test_LS_G1e-4_000/test_LS_G1e-4_000.o27502969"  , tag=0, caption="23"),\
    optresults(fi_name="test_LS_G1e-4_001/test_LS_G1e-4_001.o27502970"  , tag=0, caption="24"),\
    optresults(fi_name="test_LS_G1e-4_002/test_LS_G1e-4_002.o27502971"  , tag=0, caption="25"),\
    optresults(fi_name="test_LS_G1e-4_003/test_LS_G1e-4_003.o27502972"  , tag=0, caption="26"),\
    optresults(fi_name="test_LS_G1e-4_004/test_LS_G1e-4_004.o27502973"  , tag=0, caption="27"),\
    optresults(fi_name="test_LS_G1e-4_005/test_LS_G1e-4_005.o27502974"  , tag=0, caption="28"),\
    optresults(fi_name="test_LS_G1e-4_006/test_LS_G1e-4_006.o27502975"  , tag=0, caption="29"),\
    optresults(fi_name="test_LS_G1e-4_007/test_LS_G1e-4_007.o27502976"  , tag=0, caption="30"),\
    optresults(fi_name="test_LS_G1e-4_008/test_LS_G1e-4_008.o27502977"  , tag=0, caption="31"),\
    optresults(fi_name="test_LS_G1e-4_009/test_LS_G1e-4_009.o27502978"  , tag=0, caption="32"),\
    optresults(fi_name="test_LS_G1e-4_010/test_LS_G1e-4_010.o27502979"  , tag=0, caption="33"),\
    optresults(fi_name="test_LS_G1e-4_011/test_LS_G1e-4_011.o27502980"  , tag=0, caption="34"),\
    optresults(fi_name="test_TR_001/test_TR_001.o27503498"              , tag=1, caption="35"),\
    optresults(fi_name="test_TR_002/test_TR_002.o27503499"              , tag=1, caption="36"),\
    optresults(fi_name="test_TR_003/test_TR_003.o27503500"              , tag=1, caption="37"),\
    #optresults(fi_name="test_TR_004/test_TR_004.o27503502"              , tag=1, caption="38"),\
    #optresults(fi_name="test_TR_005/test_TR_001.o27503503"              , tag=1, caption="39"),\
                            ]

plt.close(2)
fig = plt.figure(num=2, figsize=(25.6 , 12.91))
ax = fig.add_subplot(111)
ax.set_position([0.05, 0.05, 0.85, 0.85])
k= 0 
for i, res in enumerate(all_results):
    to_plot = res.misfit()
    try:
        #ax.plot(range(len(to_plot[0][:,0])), to_plot[0][:,0]/to_plot[0][0,0],  marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_results))(i), label=str(res.caption)) 
                                                                                                                                                                                                          ax.plot(to_plot[0][:,1], to_plot[0][:,0]/to_plot[0][0,0],  marker='o', markersize=4, picker=10,color=cm.get_cmap(name='gist_rainbow', lut=len(all_results))(i), label=res.caption) 
    except:
        print(res.fi_name)
        pass
#ax.set_ylim(1e-6,1.1)
ax.set_yscale('log')
ax.legend(loc='best')
ax.grid()
fig.canvas.callbacks.connect('pick_event', on_pick)
fig.show()
