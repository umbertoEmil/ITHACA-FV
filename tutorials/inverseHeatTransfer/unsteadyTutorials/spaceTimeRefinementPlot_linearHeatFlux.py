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

folderLinear = "./linearBasisTest/resultsForPlots/linearTrueHeat/"

deltaX = np.loadtxt(folderLinear + "deltaX.txt") 
deltaX = np.prod(deltaX, axis=1)
deltaX = np.flip(deltaX)

deltaT = np.loadtxt(folderLinear + "deltaT.txt")

heatFluxRelErr_L2_linear = np.loadtxt(folderLinear + "heatFluxRelErr_L2_list.txt")
heatFluxRelErr_Linf_linear = np.loadtxt(folderLinear + "heatFluxRelErr_Linf_list.txt")

#heatFluxRelErr_L2 han on each row the values for a fixed mesh
heatFluxRelErr_L2_linear = np.reshape(heatFluxRelErr_L2_linear, (deltaX.size, deltaT.size))
heatFluxRelErr_L2_linear = np.flip(heatFluxRelErr_L2_linear, 0)

heatFluxRelErr_Linf_linear = np.reshape(heatFluxRelErr_Linf_linear, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_linear = np.flip(heatFluxRelErr_Linf_linear, 0)

folderConstant = "./constantBasisTest/resultsForPlots/linearTrueHeat/"

heatFluxRelErr_L2_constant = np.loadtxt(folderConstant + "heatFluxRelErr_L2_list.txt")
heatFluxRelErr_Linf_constant = np.loadtxt(folderConstant + "heatFluxRelErr_Linf_list.txt")

heatFluxRelErr_L2_constant = np.reshape(heatFluxRelErr_L2_constant, (deltaX.size, deltaT.size))
heatFluxRelErr_L2_constant = np.flip(heatFluxRelErr_L2_constant, 0)

heatFluxRelErr_Linf_constant = np.reshape(heatFluxRelErr_Linf_constant, (deltaX.size, deltaT.size))
heatFluxRelErr_Linf_constant = np.flip(heatFluxRelErr_Linf_constant, 0)


#############################################################
fig, axes =  plt.subplots(1)
i=0
j = len(deltaX)
for dx in deltaX:
    line, = plt.loglog(deltaT, heatFluxRelErr_L2_linear[i,:],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_L2_linear[i,:],'-', color = line.get_color(), label = "Mesh " + str(j), linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_L2_constant[i,:],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel(r"$\Delta t$ $[s]$")
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$")
plt.xlim(deltaT[0],deltaT[-1])

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["s", "o"]]
plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], fontsize=25, ncol = 2, bbox_to_anchor=(0., 1.1), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)


plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dx in deltaX:
    line, = plt.loglog(deltaT, heatFluxRelErr_Linf_linear[i,:],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_Linf_linear[i,:],'-', color = line.get_color(), label = "Mesh " + str(j), linewidth = 2, markersize = 15)
    plt.loglog(deltaT, heatFluxRelErr_Linf_constant[i,:],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel(r"$\Delta t$ $[s]$")
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$")

plt.xlim(deltaT[0],deltaT[-1])

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["s", "o"]]
plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], fontsize=25, ncol = 2, bbox_to_anchor=(0., 1.1), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dt in deltaT:
    line, = plt.loglog(deltaX, heatFluxRelErr_L2_linear[:,i],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_L2_linear[:,i],'-', color = line.get_color(), label = r"$\Delta t$ = " + str(dt), linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_L2_constant[:,i],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel("Mesh size")
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$")

plt.xlim(deltaX[0],deltaX[-1])

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["s", "o"]]
plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], fontsize=25, ncol = 2, bbox_to_anchor=(0., 1.1), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)

##############################################################
fig, axes =  plt.subplots()
i=0
j = len(deltaX)
for dt in deltaT:
    line, = plt.loglog(deltaX, heatFluxRelErr_Linf_linear[:,i],'-o', linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_Linf_linear[:,i],'-', color = line.get_color(), label = r"$\Delta t$ = " + str(dt), linewidth = 2, markersize = 15)
    plt.loglog(deltaX, heatFluxRelErr_Linf_constant[:,i],'--s', color = line.get_color(), linewidth = 2, markersize = 15)
    i = i + 1
    j = j - 1

plt.xlabel("Mesh size")
plt.ylabel(r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", fontsize=30)

plt.xlim(deltaX[0],deltaX[-1])

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["s", "o"]]
plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], fontsize=25, ncol = 2, bbox_to_anchor=(0., 1.1), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.grid(True)


plt.show()
