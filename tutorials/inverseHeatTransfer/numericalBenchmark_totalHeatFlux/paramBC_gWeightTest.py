import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relError_paramBC = np.loadtxt("coarseMesh/ITHACAoutput/paramBC_gWeightTest/totalHeatRelError_mat.txt")
relError_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/paramBC_gWeightTest/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/paramBC_gWeightTest/relError_LinfNorm_mat.txt")
gWeights = np.loadtxt("coarseMesh/ITHACAoutput/paramBC_gWeightTest/gWeightsVector_mat.txt")

markerS = 7

f = plt.figure(3,figsize=(12,8))
plt.loglog(gWeights, relError_L2norm, 'bo', markersize=markerS, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.loglog(gWeights, relError_LinfNorm, 'kv', markersize=markerS, label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel(r'Weight of the total heat flux measurement, $p_g\ [\frac{K^2}{W^2}]$', fontsize=25)
ax1 = plt.axes()
#ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))


plt.grid()
plt.title("Parameterized BC (LU)", fontsize=25)
plt.legend(fontsize=25)
plt.show()
