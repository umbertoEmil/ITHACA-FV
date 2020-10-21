import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pylab as pylab
import numpy as np
import sys

def Dt2CFL(Dt):
    diffusivity = 0.0333333333333333
    Ncells = 40
    edgeLength = 1.0 / Ncells
    return Dt * diffusivity * 3 / (edgeLength * edgeLength)


sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relErr_L2 = np.loadtxt("./timeConvL2error.txt")
relErr_Linf = np.loadtxt("./timeConvLinfError.txt")
Ncells = 40
deltaT = np.loadtxt("./timeConvDt.txt")
edgeLength = 1.0 / Ncells
CFL = np.loadtxt("./timeConvCFL.txt")


fig = plt.figure(20,figsize=(12,8))
ax1 = plt.axes(xscale='log', yscale='log')

ax1.plot(deltaT, relErr_L2, "bo", markersize=15, label=r"$|| \epsilon ||_{L^2(\Omega \times (0, t_f])}$")
ax1.plot(deltaT, relErr_Linf, "kv", markersize=15,label=r"$|| \epsilon ||_{L^\infty(\Omega \times (0, t_f])}$") 


#plt.title(r"$Mesh = $", fontsize=25)
plt.xlabel(r'$\Delta t$ [$s$]', fontsize=25)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc='best', fontsize=25)

ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
plt.xscale('log')
ax2.set_xlabel('Courant number', fontsize=25)

# get the primary axis x tick locations in plot units
xtickloc = ax1.get_xticks() 
# set the second axis ticks to the same locations
ax2.set_xticks(xtickloc)
# calculate new values for the second axis tick labels, format them, and set them
x2labels = ['{:.3g}'.format(x) for x in Dt2CFL(xtickloc)]
ax2.set_xticklabels(x2labels)
# force the bounds to be the same
ax2.set_xlim(ax1.get_xlim())



plt.show()
