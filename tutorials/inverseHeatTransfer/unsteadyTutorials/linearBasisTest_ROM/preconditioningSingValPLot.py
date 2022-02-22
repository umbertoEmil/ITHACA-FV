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

singVal = np.loadtxt("./ITHACAoutput/regularization/preconditioning/singularValues_mat.txt")
singVal_prec = np.loadtxt("./ITHACAoutput/regularization/preconditioning/singularValues_precond_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.semilogy(singVal, "X",linewidth = 2, label = "singVal")
plt.semilogy(singVal_prec, "o",linewidth = 2, label = "singVal_prec")

plt.grid()
plt.legend()

print("Original conditioning = ")
print singVal[0] / singVal[len(singVal) - 1]

print("New conditioning = ")
print singVal_prec[0] / singVal_prec[len(singVal) - 1]

##############################################################
plt.show()

