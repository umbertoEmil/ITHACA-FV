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

relErr_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/CGtest/relError_L2norm_mat.txt")
relErr_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/CGtest/relError_LinfNorm_mat.txt")
J = np.loadtxt("coarseMesh/costFunctionFull_mat.txt")

print relErr_L2norm[-1]
print relErr_LinfNorm[-1]

iterations = np.linspace(1, len(relErr_LinfNorm), len(relErr_LinfNorm))
iterationTicks = np.linspace(1, 100, 11)

f = plt.figure(2,figsize=(12,8))
#plt.semilogy(iterations, relErr_L2norm,'b-', linewidth=2)
plt.plot(iterations, relErr_L2norm,'b-', linewidth=2, label = r'$||\epsilon_{rel}||_{L^2(\Gamma_{s_in})}$')
plt.plot(iterations, relErr_LinfNorm,'k--', linewidth=2, label = r'$||\epsilon_{rel}||_{L^\infty(\Gamma_{s_in})}$')
plt.ylim(0,1)
plt.xlim(1,len(iterations))

#plt.xticks(iterationTicks)

#plt.semilogy(xAxis, singVal, "o", markersize=15)
plt.legend(loc='best', fontsize=25)
plt.xlabel("Iteration", fontsize=25)
plt.grid()

f = plt.figure(3,figsize=(12,8))
plt.semilogy(iterations, J,'b-', linewidth=2)
plt.xlabel("Iteration", fontsize=25)
plt.ylabel(r"Cost function, $J$", fontsize=25)
plt.xlim(1,len(iterations))
plt.grid()


plt.show()
