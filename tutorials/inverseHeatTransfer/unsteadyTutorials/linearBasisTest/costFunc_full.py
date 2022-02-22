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

costFunction = np.loadtxt("./ITHACAoutput/testInverse/costFunction_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.semilogy(costFunction, "b", linewidth = 2)

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Cost function, $S^k[\mathbf{w}^k]$', fontsize=25)
#plt.xlim(0,5)
plt.grid()

##############################################################
plt.show()

