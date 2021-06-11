import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

singularValues = np.loadtxt("coarseMesh/ITHACAoutput/offlineParamBC/singularValues_mat.txt")

f = plt.figure(3,figsize=(12,8))
plt.semilogy(singularValues, 'bo', markersize=15, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.ylabel('Singular values', fontsize=25)

plt.grid()
plt.show()
