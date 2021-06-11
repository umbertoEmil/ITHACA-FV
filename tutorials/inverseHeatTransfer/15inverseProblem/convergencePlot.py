import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

iteration = list()
J= list()
gradJ = list()
L2norm = list()
LinftyNorm = list()


with open("log.dat", "r") as data:
    for line in data:
	if(line[0]!="#"):
	    floats = [float(x) for x in line.split()]
	    iteration.append(int(floats[0]))
	    J.append(floats[1])
	    gradJ.append(floats[2])
	    L2norm.append(floats[3])
	    LinftyNorm.append(floats[4])


f = plt.figure(1)
ax = plt.axes(yscale='log')
#ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
plt.plot(iteration, J, linewidth=2.0)
plt.plot(iteration, gradJ, linewidth=2.0)
plt.ylabel(r'$J\    [K^2 / m^2]$' )
plt.xlabel('Iteration')
#plt.yscale('log')
plt.xlim(iteration[0], iteration[-1])
#plt.yticks(np.arange(J[-1], J[0]))
plt.grid()
plt.savefig('J.eps', bbox_inches='tight')


g = plt.figure(2)
plt.plot(iteration, L2norm, linewidth=2, label=r"$||\epsilon||_{L^2(\Sigma_{in})}$")
plt.plot(iteration, LinftyNorm, linewidth=2, label=r"$||\epsilon ||_{L^\infty(\Sigma_{in})}$")
leg = plt.legend(loc='best')
plt.xlim(iteration[0], iteration[-1])
plt.ylim(0, 2.0)
#plt.ylabel(r'$L2$')
plt.xlabel('Iteration')
plt.grid()
plt.savefig('errorNorms.eps', bbox_inches='tight')
plt.show()
