import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker


params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

deltaX = np.loadtxt("deltaX.txt") 
deltaX = np.prod(deltaX, axis=1)
deltaT = np.loadtxt("deltaT.txt")
heatFluxRelErr_L2 = np.loadtxt("heatFluxRelErr_L2_list.txt")
heatFluxRelErr_Linf = np.loadtxt("heatFluxRelErr_Linf_list.txt")

#heatFluxRelErr_L2 han on each row the values for a fixed mesh
heatFluxRelErr_L2 = np.reshape(heatFluxRelErr_L2, (deltaX.size, deltaT.size))

heatFluxRelErr_Linf = np.reshape(heatFluxRelErr_Linf, (deltaX.size, deltaT.size))

##############################################################
fig = plt.figure(1,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_L2[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=25)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=25)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(2,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_Linf[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=25)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=25)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(3,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_L2[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=25)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=25)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(4,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_Linf[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=25)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=25)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)



plt.show()
