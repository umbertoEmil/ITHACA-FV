import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 25,
          'figure.figsize': (12, 12),
         'axes.labelsize': 30,
         'axes.titlesize': 30,
         'xtick.labelsize':25,
         'ytick.labelsize':25,
         'figure.autolayout': True}
pylab.rcParams.update(params)

deltaT = [0.1, 0.2, 0.25, 0.5]
deltaT_name = ["01", "02", "025", "05"]
mesh = 3

folders = []

for tI in deltaT_name:
    folders.append("Mesh" + str(mesh) + "_Dt_" + str(tI) + "/")
    print folders

fig, axes =  plt.subplots(2, 1, sharex=True)
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0)

#First plot
i = 0
for fI in folders:
    costFunctionParameter = np.loadtxt(fI + "costFunctionParameter.txt")
    measurementsDiscrepancy_mean = np.loadtxt(fI + "measurementsDiscrepancy_mean_list.txt")

    axes[0].loglog(costFunctionParameter, measurementsDiscrepancy_mean,'-o', markersize = 6, linewidth = 2, label=r"$\Delta t = $" + str(deltaT[i]) + " s")
    i = i + 1

import itertools
for l, ms in zip(axes[0].lines, itertools.cycle('os>^+*')):
    l.set_marker(ms)

#plt.xlim(1e-16,1e-6)
axes[0].set_ylim(1e-1,1e15)
#plt.legend()
#plt.title("Mesh " + mesh)
axes[0].grid(True)
#leg = plt.legend(loc='best')
#plt.xlabel(r"Cost function parameter, $p_g \left[ \frac{K^2}{W^2} \right]$")
axes[0].set_ylabel(r"$mean_k(|| \mathbf{T}_s^k - \hat{\mathbf{T}}^k||^2)$")

#Second plot
i = 0
for fI in folders:
    costFunctionParameter = np.loadtxt(fI + "costFunctionParameter.txt")
    heatFluxL2norm_mean = np.loadtxt(fI + "heatFluxL2norm_mean_list.txt")

    axes[1].loglog(costFunctionParameter, heatFluxL2norm_mean,'-o', markersize = 6, linewidth = 2, label=r"$\Delta t =$" + str(deltaT[i]) + " s")
    i = i + 1

import itertools
for l, ms in zip(axes[1].lines, itertools.cycle('os>^+*')):
    l.set_marker(ms)

#fig, axs = plt.subplots(3, 1, sharex=True)
## Remove horizontal space between axes
#fig.subplots_adjust(hspace=0)
#
## Plot each graph, and manually set the y tick values
#axs[0].plot(t, s1)
#axs[0].set_yticks(np.arange(-0.9, 1.0, 0.4))
#axs[0].set_ylim(-1, 1)
#
#axs[1].plot(t, s2)
#axs[1].set_yticks(np.arange(0.1, 1.0, 0.2))
#axs[1].set_ylim(0, 1)


plt.xlim(1e-16,1e-6)
plt.ylim(1e4,1e22)
plt.legend()
#plt.title("Mesh " + mesh)
plt.grid(True)
#leg = plt.legend(loc='best')
plt.xlabel(r"Cost function parameter, $p_g \left[ \frac{K^2}{W^2} \right]$")
plt.ylabel(r"$mean_k(||g^k||_{L^2(\Gamma_{s_{in}})})$")

plt.show()

