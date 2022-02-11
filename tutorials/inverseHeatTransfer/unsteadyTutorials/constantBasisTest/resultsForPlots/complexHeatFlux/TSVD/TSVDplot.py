import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker


params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

alpha = np.loadtxt("TSVDalpha.txt")
heatFluxRelErr_L2 = np.loadtxt("heatFluxRelErr_L2_TSVD.txt")
heatFluxRelErr_Linf = np.loadtxt("heatFluxRelErr_Linf_TSVD.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.semilogy(alpha, heatFluxRelErr_L2,'-o', label = "L2", markersize = 15)
plt.xlabel(r"$\alpha_{TSVD}$", fontsize=25)

plt.legend(loc='upper right')
plt.grid(True)

plt.show()
