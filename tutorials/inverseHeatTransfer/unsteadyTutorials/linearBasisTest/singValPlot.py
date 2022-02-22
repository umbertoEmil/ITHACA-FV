import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
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

singVal = np.loadtxt("./ITHACAoutput/testInverse/singularValues_mat.txt")

##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.semilogy(singVal, "o",linewidth = 2, label = r"$\sigma_i$")

plt.ylabel("Singular Values", fontsize=25)
plt.grid()

##############################################################
plt.show()

