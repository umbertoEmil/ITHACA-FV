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

forderLU = "./fullPivLU/"
noiseLevel = np.loadtxt(forderLU + "noiseLevel.txt")

heatFluxRelErr_L2 = np.loadtxt(forderLU + "heatFluxRelErr_L2_max_list.txt")

Ntests = heatFluxRelErr_L2.size / noiseLevel.size

print(Ntests)

relErr_L2norm =         np.empty([len(noiseLevel)])
relErr_L2normMin =      np.empty([len(noiseLevel)])
relErr_L2normMax =      np.empty([len(noiseLevel)])

for i in range(len(noiseLevel)):
    vec = heatFluxRelErr_L2[i * Ntests:i * Ntests + Ntests - 1]
    relErr_L2norm[i] = vec.mean()
    relErr_L2normMin[i] = np.quantile(vec, 0.1)
    relErr_L2normMax[i] = np.quantile(vec, 0.9)

###############################################################
fig = plt.figure(1)
plt.yscale('log')
plt.xscale('log')
plt.errorbar(noiseLevel, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normMax - relErr_L2norm], markersize=15,fmt='bo-', capsize=12, label = 'LU')

################################################################
#forderTSVD = "./TSVD/regPamameter10/"
#noiseLevel = np.loadtxt(forderTSVD + "noiseLevel.txt")
#
#heatFluxRelErr_L2 = np.loadtxt(forderTSVD + "heatFluxRelErr_L2_mean_list.txt")
#
#Ntests = heatFluxRelErr_L2.size / noiseLevel.size
#
#relErr_L2norm =         np.empty([len(noiseLevel)])
#relErr_L2normMin =      np.empty([len(noiseLevel)])
#relErr_L2normmean =      np.empty([len(noiseLevel)])
#
#for i in range(len(noiseLevel)):
#    vec = heatFluxRelErr_L2[i * Ntests:i * Ntests + Ntests - 1]
#    relErr_L2norm[i] = vec.mean()
#    relErr_L2normMin[i] = np.quantile(vec, 0.1)
#    relErr_L2normmean[i] = np.quantile(vec, 0.9)
#plt.errorbar(noiseLevel, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normmean - relErr_L2norm], markersize=15,fmt='ks--', capsize=12, label = r'TSVD, $\alpha_{TSVD} = 10$')
#
################################################################
#forderTSVD = "./TSVD/regPamameter8/"
#noiseLevel = np.loadtxt(forderTSVD + "noiseLevel.txt")
#
#heatFluxRelErr_L2 = np.loadtxt(forderTSVD + "heatFluxRelErr_L2_mean_list.txt")
#
#Ntests = heatFluxRelErr_L2.size / noiseLevel.size
#
#relErr_L2norm =         np.empty([len(noiseLevel)])
#relErr_L2normMin =      np.empty([len(noiseLevel)])
#relErr_L2normmean =      np.empty([len(noiseLevel)])
#
#for i in range(len(noiseLevel)):
#    vec = heatFluxRelErr_L2[i * Ntests:i * Ntests + Ntests - 1]
#    relErr_L2norm[i] = vec.mean()
#    relErr_L2normMin[i] = np.quantile(vec, 0.1)
#    relErr_L2normmean[i] = np.quantile(vec, 0.9)
#plt.errorbar(noiseLevel, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normmean - relErr_L2norm], markersize=15,fmt='g^-.', capsize=12, label = r'TSVD, $\alpha_{TSVD} = 8$')
#
#
#
#plt.xlabel('Noise standard dev.')
#plt.ylabel(r'$mean_{t\in(0,t_f]} ||e_{rel}||_{L^2(\Gamma_{s_{in}})} $')
##plt.title('LU w. full pivoting', fontsize=25)
#
#plt.grid(True, which="both", ls="-")
#plt.legend()
#
#plt.ylim(1e-3,1e0)
#plt.xlim(noiseLevel[0],noiseLevel[-1])
plt.show()

