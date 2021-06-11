import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import sys
sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relErr_L2 = np.load("./relErrL2.npy")
relErr_Linf = np.load("./relErrLinf.npy")
Ncells = np.load("./Ncells.npy")
deltaT = np.load("./deltaT.npy")




fig = plt.figure(1,figsize=(12,8))
for column in relErr_L2:
    plt.semilogy(Ncells, column, linewidth = 2)

plt.xlabel('Cells per axis', fontsize=25)
plt.grid()
plt.legend()

plt.show()
