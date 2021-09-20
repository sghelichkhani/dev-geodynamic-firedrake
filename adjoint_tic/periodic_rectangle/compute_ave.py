import numpy as np
import h5py
import matplotlib.pyplot as plt


hf = h5py.File('FinalState.h5', 'r')
temperature = np.array(hf.get('fields/Temperature'))
hf.close()


hf = h5py.File('FinalStateYaxis.h5', 'r')
y_coords = np.array(hf.get('fields/Temperature'))
hf.close()


res = np.polynomial.legendre.Legendre.fit(y_coords,temperature, 35) 


y_vals = np.linspace(0, 1, 200)

plt.close(1)
fig = plt.figure(num=1)
ax = fig.add_subplot(111)
ax.scatter(y_coords, temperature)
ax.plot(y_vals, res(y_vals), color='k')
fig.show()

with open('radial.out', mode='w') as fi:
    for i in np.linspace(0, 1, 160):
        fi.write(f"{i}, {res(i)}\n")
    
