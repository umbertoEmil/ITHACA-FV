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

AxbNorm = np.loadtxt("./ITHACAoutput/regularization/conjugateGradient/Lcurve/AxbNorm_mat.txt")
xNorm = np.loadtxt("./ITHACAoutput/regularization/conjugateGradient/Lcurve/xNorm_mat.txt")
CGsteps = np.loadtxt("./ITHACAoutput/regularization/conjugateGradient/Lcurve/CGsteps_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.loglog(AxbNorm,xNorm, "bo-", linewidth = 2)

plt.ylabel('||x||', fontsize=25)
plt.xlabel('||Ax-b||', fontsize=25)
plt.grid()


##############################################################
plt.show()

