import numpy as np
import matplotlib.pyplot as plt


# Set up geometry:
rmin, rmax      = 1.22, 2.22
ncells, nlayers = 256, 50


# an aray of radius
rshl = np.linspace(rmin, rmax, nlayers)

# initiating layer heights with 1.
resolution_func = np.ones((nlayers))

# doubling resolution in the upper mantle above 700 km
#resolution_func[rshl>= (rmax -1000/6370)] *= 0.5

# A gaussian shaped function 
def gaussian(center, c):
    return res_amplifier*np.exp(-(rshl-center)**2/(2*c**2))

for idx, r_0 in enumerate([rmin, rmax, rmax - 660/6370]):
    # gaussian radius for the gassuan function
    c= 0.20
    # How different is the high res area from low res
    res_amplifier = 5.
    if idx ==2:
        res_amplifier = 2.
        c = 0.005
    resolution_func *=  1/(1+gaussian(center=r_0, c=c))

resolution_func = (rmax-rmin) * resolution_func/np.sum(resolution_func)

plt.close(1)
fig = plt.figure(num=1)
ax = fig.add_subplot(111)
ax.plot(resolution_func)
fig.show()


