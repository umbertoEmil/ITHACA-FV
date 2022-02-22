import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
sys.path.insert(0, "./")

#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

folder = "./ITHACAoutput/testInverse/"

t =                         np.loadtxt(folder + "timeVec_mat.txt")
probeCenter_gRec =          np.loadtxt(folder + "probeHeat_rec_mat.txt")
probeCenter_trueHeatFlux =  np.loadtxt(folder + "probeHeat_true_mat.txt")

heatFluxRelErr_L2_FOM = np.loadtxt(folder + "heatFluxRelErr_L2_mat.txt")
heatFluxRelErr_Linf_FOM = np.loadtxt(folder + "heatFluxRelErr_Linf_mat.txt")

##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.plot(t, probeCenter_trueHeatFlux, "r", linewidth = 2, label = r"$g_t$")
plt.plot(t, probeCenter_gRec, "bo-", linewidth = 2, label = r"$g_{FOM}$")


plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux $[W/m^2]$', fontsize=25)
plt.xlim(0,50)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig, axes =  plt.subplots(figsize=(12,8))
plt.semilogy(t,heatFluxRelErr_L2_FOM,   "bo-", linewidth = 2,markersize = 6, label = r"$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$")
plt.semilogy(t,heatFluxRelErr_Linf_FOM, "ks--", linewidth = 2,markersize = 6, label = r"$||e_{rel}||_{L^\infty(\Gamma_{s_{in}})}$")

plt.xlabel('Time [s]', fontsize=25)
plt.xlim(0,50)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()

leg = plt.legend(loc='best')
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
#plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], ncol = 2, bbox_to_anchor=(0.001, 0.06),loc=2, borderaxespad=0.)
plt.legend(handles=h, labels=[r"$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$", r"$||e_{rel}||_{L^\infty(\Gamma_{s_{in}})}$"], ncol = 2, bbox_to_anchor=(0., 1.056), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)


##############################################################
plt.show()

