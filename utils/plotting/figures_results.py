import matplotlib
import numpy as np
import glob
import os
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm
import imageio
import sys 


def show_image(fi, ax):
    image =  imageio.imread(fi)
    ax.imshow(image[200:-40,50:-50,:])

y_0   = 0.72
x_0   = 0.10
del_x = 0.135/1.5
del_y = 0.3/1.5
mar_y = 0.01
mar_x = 0.00
num_hor = 10

# what should show up as the type of simulation
caption='ScipyInterpolation'
# This is the path that all the figures are stored in
path_simu = '/Users/sghelichkhani/Workplace/develope_G-ADOPT/utils/plotting/PNGs/rectangle_extruded_ns_HighRa_Scipy/BFGS_linesearchCubic/'
# This is the path with all the initial/final figures
fwdref_path = '/Users/sghelichkhani/Workplace/develope_G-ADOPT/utils/plotting/PNGs/rectangle_extruded_ns_HighRa_Scipy/FWD_refmodel/'
# Number of iterations to jump
iter_jump = 60
# What is the index of the reference initial and final condition
init_index = 0 
fin_index = int(np.array([float(i[i.rfind('_')+1:i.rfind('.png')]) for i in glob.glob(os.path.join(fwdref_path, "*png"))]).max())
# this is the number of maximum figures out ther that can be read in
max_iter = int(np.array([float(i[i.rfind('_')+1:i.rfind('.png')]) for i in glob.glob(os.path.join(path_simu, "*png"))]).max())
max_iter= 241

plt.close(1)
fig = plt.figure(num=1, figsize=(24, 12))
# ax will be the reference initial condition
ax = fig.add_subplot(111)
ax.set_position([x_0-0.2*del_x, y_0, del_x, del_y])


ax.text(0.50, 1.08, 'What we want to retrieve:\nReference Initial Condition', rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=ax.transAxes,\
        bbox={'facecolor':(1.0, 1.0, 1.0), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})

# Plotting the initial reference figure 
show_image(os.path.join(fwdref_path, str('fig_%3.3i.png' %init_index)), ax)
ax.set_axis_off()

# i will be used to position the ax vertically and horizontally
i = 0
# bx is a list of all the figures
bx = []
for iteration in range(0, max_iter, iter_jump):
    # position the ax based on "i"
    axis_position = [x_0+(np.mod(i+1, num_hor-1))*(del_x+ mar_x), y_0-int((i+1)/(num_hor-1))*(del_y+mar_y), del_x, del_y]
    # add the ax to the list
    bx.append(fig.add_subplot(111))
    bx[-1].set_position(axis_position)
    # show the image
    show_image(os.path.join(path_simu, str('fig_init_%2.2i.png' %iteration)), bx[-1])
    # Text for the figure (Which iteration is this one?) 
    bx[-1].text(0.50, 1.05, str('Iteration: %i' %iteration), rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=bx[-1].transAxes,\
            bbox={'facecolor':(0.9, 0.9, 0.9), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})
    bx[-1].set_axis_off()
    i+=1

# Now we go to the second half of the figure (FWD conditions: reconstructed and reference)
y_0 -= 1.1*del_y+mar_y 
cx = fig.add_subplot(111)
cx.set_position([x_0-0.2*del_x, y_0, del_x, del_y])

# Reference Final figure
cx.text(0.50, 1.08, 'What we have:\nReference Final Condition', rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=cx.transAxes,\
        bbox={'facecolor':(1.0, 1.0, 1.0), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})

show_image(os.path.join(fwdref_path, str('fig_%3.3i.png' %fin_index)), cx)
cx.set_axis_off()

# These are the reconstructions of the final conditions 
i = 0
for iteration in range(0, max_iter, iter_jump):
    axis_position = [x_0+(np.mod(i+1, num_hor-1))*(del_x+ mar_x), y_0-int((i+1)/(num_hor-1))*(del_y+mar_y), del_x, del_y]
    bx.append(fig.add_subplot(111))
    bx[-1].set_position(axis_position)
    show_image(os.path.join(path_simu, str('fig_fin_%2.2i.png' %iteration)), bx[-1])
    #bx[-1].text(0.50, 1.05, str('Iteration: %i' %iteration), rotation=0, style='italic', ha='center', va='center', fontsize=12, transform=bx[-1].transAxes,\
    #        bbox={'facecolor':(0.9, 0.9, 0.9), 'alpha':1.0, 'lw':2.0,'pad':4,'edgecolor':'k'})
    bx[-1].set_axis_off()
    i+=1

cx.annotate('', xy=(1.12,  0.0), xycoords='axes fraction', xytext=(1.12, +2.40),\
                           arrowprops=dict(arrowstyle="-", ls='dashed', color='k', linewidth=1, mutation_scale=20))
cx.annotate('', xy=(1.1,  0.0), xycoords='axes fraction', xytext=(1.1, +2.40),\
                           arrowprops=dict(arrowstyle="-", ls='dashed', color='k', linewidth=1, mutation_scale=20))


fig.show()
fig_name="RhodriAuscope.png"#str('Tests_28July/%s.png'  %caption)
print(fig_name)
fig.savefig(fig_name, dpi=100, bbox_inches='tight', pad_inches=0.1)
