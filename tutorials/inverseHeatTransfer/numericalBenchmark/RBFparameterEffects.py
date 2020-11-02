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

singVal = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCtest_RBFparameter/singularValues_mat.txt")
relErr_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCtest_RBFparameter/relError_L2norm_mat.txt")
relErr_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCtest_RBFparameter/relError_LinfNorm_mat.txt")
condNumber = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCtest_RBFparameter/condNumber_mat.txt")

singVal = singVal / singVal[0,:][None,:]

xAxis = np.linspace(1, singVal.shape[0], num = singVal.shape[0]) 

shapePar = np.array([100, 10, 1, 0.3, 0.1, 0.03, 0.01, 0.0033, 0.001])
labels=[r'$\eta = 100$', r'$\eta = 10$', r'$\eta = 1$', r'$\eta = 0.3$', r'$\eta = 0.1$', r'$\eta = 0.03$', r'$\eta = 0.01$', r'$\eta = 0.0033$', r'$\eta = 0.001$']




f = plt.figure(2,figsize=(14,8))
for i in range(singVal.shape[1]):
    plt.semilogy(xAxis,singVal[:,i],'o', markersize=15, label=labels[i])



#plt.semilogy(xAxis, singVal, "o", markersize=15)
plt.legend(loc='best', fontsize=15)
plt.ylabel("Normalized singular values", fontsize=25)
plt.grid()



fig, ax1 = plt.subplots(figsize=(14,8))

color = 'tab:red'
ax1.set_xlabel(r'RBF shape parameter, $\eta$', fontsize=25)
ax1.set_ylabel('Relative error', fontsize=25, color=color)
ax1.loglog(shapePar, relErr_L2norm, 'o', markersize=15, label = r'$L^2-norm$', color=color)
ax1.loglog(shapePar, relErr_LinfNorm, 'v', markersize=15, label = r'$L^\infty-norm$', color=color)
ax1.tick_params(axis='y', labelcolor=color)

plt.grid()
plt.legend(loc=3, fontsize=25)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Condition number',fontsize=25, color=color)  # we already handled the x-label with ax1
ax2.loglog(shapePar, condNumber, 'X', markersize=15, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped



#plt.loglog(shapePar, relErr_L2norm, 'bo', markersize=15, label = r'$L^2-norm$')
#plt.loglog(shapePar, relErr_LinfNorm, 'kv', markersize=15, label = r'$L^\infty-norm$')
#plt.loglog(shapePar, condNumber, 'gX', markersize=15, label = r'$L^\infty-norm$')
#plt.xlabel(r'RBF shape parameter, $\rho$', fontsize=25)
#plt.ylabel('Relative error', fontsize=25)
#plt.grid()
#plt.legend(loc='best', fontsize=25)

plt.show()