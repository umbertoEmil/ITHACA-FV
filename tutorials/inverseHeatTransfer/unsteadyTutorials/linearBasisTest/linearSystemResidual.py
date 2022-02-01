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

heatFluxRelErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/heatFluxRelErr_L2_mat.txt")
heatFluxRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/heatFluxRelErr_Linf_mat.txt")
projHeatFluxRelErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/projHeatFluxRelErr_L2_mat.txt")
projHeatFluxRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/projHeatFluxRelErr_Linf_mat.txt")
projRecHeatFluxRelErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/projRecHeatFluxRelErr_L2_mat.txt")
projRecHeatFluxRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/projRecHeatFluxRelErr_Linf_mat.txt")
timeVec = np.loadtxt("./ITHACAoutput/testInverse/timeVec_mat.txt")
probeHeat_rec = np.loadtxt("./ITHACAoutput/testInverse/probeHeat_rec_mat.txt")
probeHeat_proj = np.loadtxt("./ITHACAoutput/testInverse/probeHeat_proj_mat.txt")
probeHeat_true = np.loadtxt("./ITHACAoutput/testInverse/probeHeat_true_mat.txt")


##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.semilogy(timeVec, projHeatFluxRelErr_L2, "o-", linewidth = 2, label = "L2")
plt.semilogy(timeVec, projHeatFluxRelErr_Linf, "o-", linewidth = 2, label = "Linf")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Projected heat flux relative error norms [$W/m$]', fontsize=25)
#plt.xlim(0,5)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(3,figsize=(12,8))
plt.semilogy(timeVec,heatFluxRelErr_L2, "o-", linewidth = 2, label = "L2")
plt.semilogy(timeVec,heatFluxRelErr_Linf, "o-", linewidth = 2, label = "Linf")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux relative error norms [$W/m$]', fontsize=25)
#plt.xlim(0,5)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(30,figsize=(12,8))
plt.semilogy(timeVec, projRecHeatFluxRelErr_L2, "o-", linewidth = 2, label = "L2")
plt.semilogy(timeVec, projRecHeatFluxRelErr_Linf, "o-", linewidth = 2, label = "Linf")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Proj heat flux relative error norms [$W/m$]', fontsize=25)
#plt.xlim(0,5h)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(4,figsize=(12,8))
plt.plot(timeVec, probeHeat_true, "b", linewidth = 2, label = r"$g_{true}$")
plt.plot(timeVec, probeHeat_rec, "k--", linewidth = 2, label = r"$g$")
plt.plot(timeVec, probeHeat_proj, "g.", linewidth = 2, label = r"$g_{proj}$")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux [$W/m^2$]', fontsize=25)
#plt.xlim(0,5)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(5,figsize=(12,8))
plt.plot(timeVec, probeHeat_true, "b", linewidth = 2, label = r"$g_t$")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux [$W/m^2$]', fontsize=25)
#plt.xlim(0,5)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()


##############################################################
plt.show()

