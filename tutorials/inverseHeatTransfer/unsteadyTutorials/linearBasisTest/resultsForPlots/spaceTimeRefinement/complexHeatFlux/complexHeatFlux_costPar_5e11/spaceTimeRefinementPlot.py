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
heatFluxRelErr_L2_max = np.loadtxt("heatFluxRelErr_L2_max_list.txt")
heatFluxRelErr_Linf_max = np.loadtxt("heatFluxRelErr_Linf_max_list.txt")

#heatFluxRelErr han on each row the values for a fixed mesh
heatFluxRelErr_L2_max = np.reshape(heatFluxRelErr_L2_max, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_max = np.reshape(heatFluxRelErr_Linf_max, (deltaX.size, deltaT.size))

heatFluxRelErr_L2_mean = np.loadtxt("heatFluxRelErr_L2_mean_list.txt")
heatFluxRelErr_Linf_mean = np.loadtxt("heatFluxRelErr_Linf_mean_list.txt")

#heatFluxRelErr han on each row the values for a fixed mesh
heatFluxRelErr_L2_mean = np.reshape(heatFluxRelErr_L2_mean, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_mean = np.reshape(heatFluxRelErr_Linf_mean, (deltaX.size, deltaT.size))

##############################################################
fig = plt.figure(1,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_L2_max[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=30)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(2,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_Linf_max[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=30)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(3,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_L2_max[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=30)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(4,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_Linf_max[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=30)
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(11,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_L2_mean[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=30)
plt.ylabel(r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(12,figsize=(12,8))
i=0
for dx in deltaX:
    plt.loglog(deltaT, heatFluxRelErr_Linf_mean[i,:],'-o', label = "mesh size = {:.1e}".format(dx), markersize = 15)
    i = i + 1
plt.xlabel(r"$\Delta t$", fontsize=30)
plt.ylabel(r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaT[0],deltaT[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(13,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_L2_mean[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=30)
plt.ylabel(r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)

##############################################################
fig = plt.figure(14,figsize=(12,8))
i=0
for dt in deltaT:
    plt.loglog(deltaX, heatFluxRelErr_Linf_mean[:,i],'-o', label = r"$\Delta t$ = {:.1e}".format(dt), markersize = 15)
    i = i + 1
plt.xlabel("Meshsize", fontsize=30)
plt.ylabel(r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=30)
plt.xlim(deltaX[0],deltaX[-1])

plt.legend(loc='upper right')
plt.grid(True)

plt.show()
