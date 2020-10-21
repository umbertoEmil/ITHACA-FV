import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
sys.path.insert(0, "./")
relErrL2norm = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")  
relErrLinfNorm = np.loadtxt("./ITHACAoutput/testInverse/relErrLinfNorm_mat.txt")  
timeSteps = np.loadtxt("./ITHACAoutput/testInverse/timeSteps_mat.txt")  
samplingTime = np.loadtxt("./ITHACAoutput/testInverse/samplingTime_mat.txt")  

#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

f = plt.figure(1,figsize=(8,6))
plt.plot(timeSteps,relErrL2norm , "b", linewidth=3.0, label=r"$||\epsilon_{rel}||_{L^2(\Gamma_{in})}$")
plt.plot(timeSteps,relErrLinfNorm , "k", linewidth=3.0, label=r"$||\epsilon_{rel} ||_{L^\infty(\Gamma_{in})}$")
plt.vlines(samplingTime, 0, 1, "r")
plt.xlabel('Time [s]', fontsize=25)
#plt.xticks(iteration)
#plt.ylim(1e-7, 1)
#plt.xticks(np.rint(np.linspace(iteration[0], iteration[-1], num=8)))
leg = plt.legend(loc='best', fontsize=25)
plt.grid()
#plt.savefig('J.eps', bbox_inches='tight')

#g = plt.figure(2,figsize=(8,6))
#plt.plot(iteration, L2norm, "b", linewidth=3, label=r"$||\epsilon_{rel}||_{L^2(\Gamma_{in})}$")
#plt.plot(iteration, LinftyNorm, "k--",linewidth=3, label=r"$||\epsilon_{rel} ||_{L^\infty(\Gamma_{in})}$")
#leg = plt.legend(loc='best', fontsize=25)
#plt.xlim(iteration[0], iteration[-1])
#plt.ylim(0, 1)
##plt.ylabel(r'$L2$')
##plt.xticks(np.rint(np.linspace(iteration[0], iteration[-1], num=8)))
#plt.xlabel('Iteration', fontsize=25)
#plt.grid()
#plt.savefig('errorNorms.eps', bbox_inches='tight')
plt.show()
