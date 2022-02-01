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

heatFluxRelErr_L2 = np.loadtxt("./ITHACAoutput/reconstructionTest/heatFluxRelErr_L2_mat.txt")
heatFluxRelErr_Linf = np.loadtxt("./ITHACAoutput/reconstructionTest/heatFluxRelErr_Linf_mat.txt")
TsolRelErr_L2 = np.loadtxt("./ITHACAoutput/reconstructionTest/TsolRelErr_L2_mat.txt")
TsolRelErr_Linf = np.loadtxt("./ITHACAoutput/reconstructionTest/TsolRelErr_Linf_mat.txt")
TrecRelErr_L2 = np.loadtxt("./ITHACAoutput/reconstructionTest/TrecRelErr_L2_mat.txt")
TrecRelErr_Linf = np.loadtxt("./ITHACAoutput/reconstructionTest/TrecRelErr_Linf_mat.txt")
TrecTsolRelErr_L2 = np.loadtxt("./ITHACAoutput/reconstructionTest/TrecTsolRelErr_L2_mat.txt")
TrecTsolRelErr_Linf = np.loadtxt("./ITHACAoutput/reconstructionTest/TrecTsolRelErr_Linf_mat.txt")
timeVec = np.loadtxt("./ITHACAoutput/reconstructionTest/timeVec_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.semilogy(timeVec, TsolRelErr_L2, "o-", linewidth = 2, label = "Solved L2")
plt.semilogy(timeVec, TsolRelErr_Linf, "o-", linewidth = 2, label = "Solved Linf")
plt.semilogy(timeVec, TrecRelErr_L2, "o-", linewidth = 2, label = "Reconstructed L2")
plt.semilogy(timeVec, TrecRelErr_Linf, "o-", linewidth = 2, label = "Reconstructed Linf")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Temperature relative error norms [$W/m$]', fontsize=25)
#plt.xlim(0,5)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.semilogy(timeVec, TrecTsolRelErr_L2, "o-", linewidth = 2, label = "L2")
plt.semilogy(timeVec, TrecTsolRelErr_Linf, "o-", linewidth = 2, label = "Linf")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Temperature relative error norms [$W/m$]', fontsize=25)
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
plt.show()
