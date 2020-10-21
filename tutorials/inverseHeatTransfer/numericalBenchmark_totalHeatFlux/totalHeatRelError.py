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
relError_CG = np.loadtxt("coarseMesh/ITHACAoutput/CG_gWeightTest/totalHeatRelError_mat.txt")
gWeights = np.loadtxt("coarseMesh/ITHACAoutput/paramBC_gWeightTest/gWeightsVector_mat.txt")

markerS = 7

f = plt.figure(3,figsize=(14,9))
plt.loglog(gWeights, relError_paramBC, 'gX', marker = "X", markersize=markerS, label = "parameterized BC")
plt.loglog(gWeights, relError_CG, 'mD', marker = 'D', markersize=markerS, label = "Alifanov's reg.")
plt.xlabel(r'Weight of the total heat flux measurement, $p_g\ [\frac{K^2}{W^2}]$', fontsize=25)
ax1 = plt.axes()
#ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0), fontsize=25)

plt.ylabel(r'$\frac{|\int_{\Gamma_{s_{in}}} g d\Gamma - \hat{G}|}{\hat{G}}$', fontsize=25)


plt.grid()
plt.legend(fontsize=25)
plt.show()
