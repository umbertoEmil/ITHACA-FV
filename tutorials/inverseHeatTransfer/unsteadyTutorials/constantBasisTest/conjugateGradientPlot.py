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

residualNorm =        np.loadtxt("./ITHACAoutput/regularization/conjugateGradient/residualNorm_mat.txt")
solutionNorm = np.loadtxt("./ITHACAoutput/regularization/conjugateGradient/solutionNorm_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.semilogy(residualNorm, "o", linewidth = 2, label = "||r||")
plt.semilogy(solutionNorm, "o", linewidth = 2, label = "||x||")
#
#plt.xlabel('Regularization parameter', fontsize=25)
#plt.ylabel('GCV function', fontsize=25)
plt.grid()
plt.legend()


##############################################################
plt.show()

