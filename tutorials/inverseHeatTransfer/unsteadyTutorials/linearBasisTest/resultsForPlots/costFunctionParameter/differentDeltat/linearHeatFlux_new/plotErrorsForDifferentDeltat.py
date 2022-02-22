import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 25,
          'figure.figsize': (12, 8),
         'axes.labelsize': 30,
         'axes.titlesize': 30,
         'xtick.labelsize':25,
         'ytick.labelsize':25,
         'figure.autolayout': True}
pylab.rcParams.update(params)

deltaT = [0.1, 0.2, 0.25, 0.5]
deltaT_name = ["01", "02", "025", "05"]
mesh = "3"

folders = []

for dt in deltaT_name:
    folders.append("Mesh" + mesh + "_Dt_" + dt + "/")
    print folders

fig, axes =  plt.subplots(1)
i = 0
for fI in folders:
    costFunctionParameter = np.loadtxt(fI + "costFunctionParameter.txt")
    heatFluxRelErr_L2_max = np.loadtxt(fI + "heatFluxRelErr_L2_max_list.txt")

    plt.loglog(costFunctionParameter, heatFluxRelErr_L2_max,'-o', markersize = 6, linewidth = 2, label=r"$\Delta t = $" + str(deltaT[i]) + " s")
    i = i + 1

import itertools
for l, ms in zip(axes.lines, itertools.cycle('os>^+*')):
    l.set_marker(ms)

plt.xlim(1e-16,1e-6)
plt.ylim(1e-4,1e14)
plt.legend()
plt.title("Mesh " + mesh)
plt.grid(True)
leg = plt.legend(loc='best')
plt.xlabel(r"Cost function parameter, $p_g \left[ \frac{K^2}{W^2} \right]$")
plt.ylabel(r"$max_{t \in (0, t_f]}\left( ||e_{rel}||_{L^2(\Gamma_{s_{in}})}\right)$")

plt.show()

