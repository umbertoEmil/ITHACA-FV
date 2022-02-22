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

costFunctionParameter = np.loadtxt("costFunctionParameter.txt")
deltaT = np.loadtxt("deltaT.txt")
deltaX = np.loadtxt("deltaX.txt")

heatFluxRelErr_L2_max = np.loadtxt("heatFluxRelErr_L2_max_list.txt")
heatFluxRelErr_Linf_max = np.loadtxt("heatFluxRelErr_Linf_max_list.txt")

heatFluxRelErr_L2_mean = np.loadtxt("heatFluxRelErr_L2_mean_list.txt")
heatFluxRelErr_Linf_mean = np.loadtxt("heatFluxRelErr_Linf_mean_list.txt")


###############################################################
fig, axes =  plt.subplots(1)

plt.loglog(costFunctionParameter, heatFluxRelErr_L2_max,'b-o', markersize = 7, linewidth = 2)
plt.loglog(costFunctionParameter, heatFluxRelErr_Linf_max,'ks', markersize = 7)
plt.loglog(costFunctionParameter, heatFluxRelErr_Linf_max,'k-', markersize = 7, linewidth = 2, label=r"$max_{t \in (0, t_f]}\left( ||e_{rel}||\right)$")

plt.loglog(costFunctionParameter, heatFluxRelErr_L2_mean,'b:o', markersize = 7, linewidth = 2)
plt.loglog(costFunctionParameter, heatFluxRelErr_Linf_mean,'ks', markersize = 7)
plt.loglog(costFunctionParameter, heatFluxRelErr_Linf_mean,'k:', markersize = 7, linewidth = 2, label=r"$mean_{t \in (0, t_f]}\left( ||e_{rel}||\right)$")
plt.xlabel(r"Cost function parameter, $p_g \left[ \frac{K^2}{W^2} \right]$", fontsize=25)

plt.legend()
#plt.title(r"$\Delta t = $" + str(deltaT) + " s , mesh = " + str(deltaX) , fontsize = 25 )
plt.grid(True)
leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
plt.legend(handles=h, labels=[r"$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$", r"$||e_{rel}||_{L^\infty(\Gamma_{s_{in}})}$"], fontsize=25, ncol = 1, bbox_to_anchor=(0.0, .25), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)

plt.xlim(1e-16,1e-6)
plt.ylim(1e-1,1e17)


plt.show()
