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

costFunctionParameter = np.loadtxt("costFunctionParameter.txt")
deltaT = np.loadtxt("deltaT.txt")
deltaX = np.loadtxt("deltaX.txt")

heatFluxL2norm_max = np.loadtxt("heatFluxL2norm_max_list.txt")
heatFluxL2norm_mean = np.loadtxt("heatFluxL2norm_mean_list.txt")


###############################################################
fig = plt.figure(1,figsize=(12,8))

plt.loglog(costFunctionParameter, heatFluxL2norm_max,'-o', label = r"max$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", markersize = 5)
plt.xlabel("Cost function parameter", fontsize=25)
plt.ylabel(r"Maximum of the relative error norms in $(0, t_f]$", fontsize=25)

plt.legend()
plt.title(r"$\Delta t = $" + str(deltaT) + " s , mesh = " + str(deltaX) , fontsize = 25 )
plt.grid(True)

###############################################################
fig = plt.figure(2,figsize=(12,8))

plt.loglog(costFunctionParameter, heatFluxL2norm_mean,'-o', label = r"mean$_{t\in (0, t_f]}\left(||e_{rel}||_{L^2(\Gamma_{s_{in}})} \right)$", markersize = 5)
plt.xlabel("Cost function parameter", fontsize=25)
plt.ylabel(r"Mean of the relative error norms in $(0, t_f]$", fontsize=25)

plt.legend()
plt.title(r"$\Delta t = $" + str(deltaT) + " s , mesh = " + str(deltaX) , fontsize = 25 )
plt.grid(True)

plt.show()
