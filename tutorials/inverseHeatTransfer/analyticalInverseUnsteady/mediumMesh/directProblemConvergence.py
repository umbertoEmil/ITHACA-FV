import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
sys.path.insert(0, "./")

directOutputFolder = "./ITHACAoutput/directProblem/"
relError_L2norm = np.loadtxt(directOutputFolder + "relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt(directOutputFolder + "relError_LinfNorm_mat.txt")
Error_L2norm = np.loadtxt(directOutputFolder + "Error_L2norm_mat.txt")
Error_LinfNorm = np.loadtxt(directOutputFolder + "Error_LinfNorm_mat.txt")
time = np.loadtxt(directOutputFolder + "timeSteps_mat.txt")

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)


f = plt.figure(1,figsize=(8,6))
ax = plt.axes(yscale='log')
#ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
plt.plot(time, Error_L2norm, "b", linewidth=3.0)
plt.plot(time, Error_LinfNorm, "k", linewidth=3.0)
plt.plot(time, relError_L2norm, "b--", linewidth=3.0)
plt.plot(time, relError_LinfNorm, "k--", linewidth=3.0)
#plt.ylabel(r'$J\    [K^2]$',fontsize=25 )
#plt.xlabel('Iteration', fontsize=25)
plt.xlim(time[0], time[-1])
#plt.xticks(iteration)
#plt.ylim(1e-7, 1)
#plt.xticks(np.rint(np.linspace(iteration[0], iteration[-1], num=8)))
plt.grid()
plt.savefig('errosPlot.eps', bbox_inches='tight')
#
#
#g = plt.figure(2,figsize=(8,6))
#plt.plot(iteration, L2norm, "b", linewidth=3, label=r"$||\epsilon_{rel}||_{L^2(\Gamma_{in})}$")
#plt.plot(iteration, LinftyNorm, "k--",linewidth=3, label=r"$||\epsilon_{rel} ||_{L^\infty(\Gamma_{in})}$")
#leg = plt.legend(loc='best', fontsize=25)
#plt.xlim(iteration[0], iteration[-1])
#plt.ylim(0, 1)
##plt.ylabel(r'$L2$')
#plt.xticks(np.rint(np.linspace(iteration[0], iteration[-1], num=8)))
#plt.xlabel('Iteration', fontsize=25)
#plt.grid()
#plt.savefig('errorNorms.eps', bbox_inches='tight')
plt.show()
