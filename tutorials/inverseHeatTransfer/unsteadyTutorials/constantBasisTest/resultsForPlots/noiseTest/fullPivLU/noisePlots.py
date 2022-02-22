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
         'ytick.labelsize':25,
         'figure.autolayout': True}
pylab.rcParams.update(params)

noiseLevel = np.loadtxt("noiseLevel.txt")

heatFluxRelErr_L2_max = np.loadtxt("heatFluxRelErr_L2_max_list.txt")
heatFluxRelErr_Linf_max = np.loadtxt("heatFluxRelErr_Linf_max_list.txt")

Ntests = heatFluxRelErr_L2_max.size / noiseLevel.size

print(Ntests)

relErr_L2norm =         np.empty([len(noiseLevel)])
relErr_L2normMin =      np.empty([len(noiseLevel)])
relErr_L2normMax =      np.empty([len(noiseLevel)])
relErr_LinfNorm =       np.empty([len(noiseLevel)])
relErr_LinfNormMin =    np.empty([len(noiseLevel)])
relErr_LinfNormMax =    np.empty([len(noiseLevel)])

for i in range(len(noiseLevel)):
    vec = heatFluxRelErr_L2_max[i * Ntests:i * Ntests + Ntests - 1]
    relErr_L2norm[i] = vec.mean()
    relErr_L2normMin[i] = np.quantile(vec, 0.1)
    relErr_L2normMax[i] = np.quantile(vec, 0.9)
    vec = heatFluxRelErr_Linf_max[i * Ntests:i * Ntests + Ntests - 1]
    relErr_LinfNorm[i] = vec.mean()
    relErr_LinfNormMin[i] = np.quantile(vec, 0.1)
    relErr_LinfNormMax[i] = np.quantile(vec, 0.9)

print relErr_L2norm

###############################################################
fig = plt.figure(1)
plt.yscale('log')
plt.xscale('log')
plt.errorbar(noiseLevel, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normMax - relErr_L2norm], markersize=15,fmt='bo-', capsize=12, label = r'$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$')
plt.errorbar(noiseLevel, relErr_LinfNorm, yerr=[relErr_LinfNorm - relErr_LinfNormMin, relErr_LinfNormMax - relErr_LinfNorm], markersize=15,fmt='kv-', capsize=12, label = r'$||e_{rel}||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel('Noise standard dev.')
plt.title(r'LU w. full pivoting')

plt.grid(True, which="both", ls="-")
plt.legend()
plt.ylim(1e-2,1e0)
plt.show()


#
#plt.loglog(noiseLevel, heatFluxRelErr_L2_max,'-o', label = r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", markersize = 15)
#plt.loglog(noiseLevel, heatFluxRelErr_Linf_max,'--s', label = r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^\infty(\Gamma_{s_{in}})} \right)$", markersize = 15)
#plt.xlabel(r"$\Delta t$", fontsize=25)
#
#plt.legend()
#plt.title("LU w. full pivoting", fontsize = 25 )
#plt.grid(True)
#
#plt.show()
