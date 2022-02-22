import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker


params = {'legend.fontsize': 25,
          'figure.figsize': (12, 8),
         'axes.labelsize': 30,
         'axes.titlesize': 30,
         'xtick.labelsize':25,
         'ytick.labelsize':25}
pylab.rcParams.update(params)

folder = "./"

deltaX = np.loadtxt(folder + "deltaX.txt") 
deltaX = np.prod(deltaX, axis=1)
deltaX = np.flip(deltaX)

deltaT = np.loadtxt(folder + "deltaT.txt")

heatFluxRelErr_L2_max = np.loadtxt(folder + "heatFluxRelErr_L2_max_list.txt")
heatFluxRelErr_Linf_max = np.loadtxt(folder + "heatFluxRelErr_Linf_max_list.txt")

#heatFluxRelErr_L2 han on each row the values for a fixed mesh
heatFluxRelErr_L2_max = np.reshape(heatFluxRelErr_L2_max, (deltaX.size, deltaT.size))
heatFluxRelErr_L2_max = np.flip(heatFluxRelErr_L2_max, 0)

heatFluxRelErr_Linf_max = np.reshape(heatFluxRelErr_Linf_max, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_max = np.flip(heatFluxRelErr_Linf_max, 0)

heatFluxRelErr_L2_mean = np.loadtxt(folder + "heatFluxRelErr_L2_mean_list.txt")
heatFluxRelErr_Linf_mean = np.loadtxt(folder + "heatFluxRelErr_Linf_mean_list.txt")

heatFluxRelErr_L2_mean = np.reshape(heatFluxRelErr_L2_mean, (deltaX.size, deltaT.size))
heatFluxRelErr_L2_mean = np.flip(heatFluxRelErr_L2_mean, 0)

heatFluxRelErr_Linf_mean = np.reshape(heatFluxRelErr_Linf_mean, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_mean = np.flip(heatFluxRelErr_Linf_mean, 0)


#############################################################
fig, axes =  plt.subplots(1)
i=0
j = len(deltaX)
for dx in deltaX:
    line, = plt.loglog(deltaT, heatFluxRelErr_L2_max[i,:],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_L2_max[i,:],'-', color = line.get_color(), label = "Mesh " + str(j), linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_L2_mean[i,:],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel(r"$\Delta t$ $[s]$")

plt.xlim(deltaT[0],deltaT[-1])
plt.ylim(1e-4, 1e17)

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
plt.legend(handles=h, labels=[r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$"], fontsize=25, ncol = 2, bbox_to_anchor=(-0.12, 1.15), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)


plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dx in deltaX:
    line, = plt.loglog(deltaT, heatFluxRelErr_Linf_max[i,:],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_Linf_max[i,:],'-', color = line.get_color(), label = "Mesh " + str(j), linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_Linf_mean[i,:],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel(r"$\Delta t$ $[s]$")

plt.xlim(deltaT[0],deltaT[-1])
plt.ylim(1e-4, 1e17)

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
plt.legend(handles=h, labels=[r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$"], fontsize=25, ncol = 2, bbox_to_anchor=(-0.12, 1.15), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dt in deltaT:
    line, = plt.loglog(deltaX, heatFluxRelErr_L2_max[:,i],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_L2_max[:,i],'-', color = line.get_color(), label = r"$\Delta t$ = " + str(dt), linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_L2_mean[:,i],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel("Mesh size")

plt.xlim(deltaX[0],deltaX[-1])
plt.ylim(1e-4, 1e17)

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
plt.legend(handles=h, labels=[r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$"], fontsize=25, ncol = 2, bbox_to_anchor=(-0.12, 1.15), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dt in deltaT:
    line, = plt.loglog(deltaX, heatFluxRelErr_Linf_max[:,i],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_Linf_max[:,i],'-', color = line.get_color(), label = r"$\Delta t$ = " + str(dt), linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_Linf_mean[:,i],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel("Mesh size")

plt.xlim(deltaX[0],deltaX[-1])
plt.ylim(1e-4, 1e17)

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
plt.legend(handles=h, labels=[r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$"], fontsize=25, ncol = 2, bbox_to_anchor=(-0.12, 1.15), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)


plt.show()
