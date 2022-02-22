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

t = np.loadtxt("./ITHACAoutput/testInverse/timeSteps_mat.txt")
heatFluxL2norm = np.loadtxt("./ITHACAoutput/testInverse/heatFluxL2norm_mat.txt")
heatFluxLinfNorm = np.loadtxt("./ITHACAoutput/testInverse/heatFluxLinfNorm_mat.txt")
redHeatFluxL2norm = np.loadtxt("./ITHACAoutput/testInverseReduced/heatFluxL2norm_mat.txt")
redHeatFluxLinfNorm = np.loadtxt("./ITHACAoutput/testInverseReduced/heatFluxLinfNorm_mat.txt")
trueHeatFluxL2norm = np.loadtxt("./ITHACAoutput/testInverse/trueHeatFluxL2norm_mat.txt")
trueHeatFluxLinfNorm = np.loadtxt("./ITHACAoutput/testInverse/trueHeatFluxLinfNorm_mat.txt")


#relErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")
#relErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/relErrLinfNorm_mat.txt")
#avgRelErr_L2 =   np.loadtxt("./ITHACAoutput/testInverse/avgRelErrL2norm_mat.txt")
#avgRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/avgRelErrLinfNorm_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.plot(t, trueHeatFluxL2norm, "b", linewidth = 2, label = r"$||g_t||_{L^2(\Gamma_{in})}$")
plt.plot(t, heatFluxL2norm, "k--", linewidth = 2, label = r"$||g||_{L^2(\Gamma_{in})}$")
plt.plot(t, redHeatFluxL2norm, "r--", linewidth = 2, label = r"$||g_{r}||_{L^2(\Gamma_{in})}$")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux $L^2-norm$ [$W$]', fontsize=25)
plt.xlim(0,50)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.plot(t, trueHeatFluxLinfNorm, "b", linewidth = 2, label = r"$||g_t||_{L^\infty(\Gamma_{in})}$")
plt.plot(t, heatFluxLinfNorm, "k--", linewidth = 2, label = r"$||g||_{L^\infty(\Gamma_{in})}$")
plt.plot(t, redHeatFluxLinfNorm, "r--", linewidth = 2, label = r"$||g_{r}||_{L^\infty(\Gamma_{in})}$")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux $L^\infty$-norm [$W$]', fontsize=25)
plt.xlim(0,50)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
plt.show()

