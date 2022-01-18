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

Ttrue = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/Ttrue_vec_mat.txt")
recT = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/recT_vec_mat.txt")
toCheckT = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/toCheckT_vec_mat.txt")
recT_relErr = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/recT_relErr_vec_mat.txt")
toCheckTrecT_relErr = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/toCheckTrecT_relErr_vec_mat.txt")
toCheckT_relErr = np.loadtxt("./ITHACAoutput/thermocoupleReconstructionTest/toCheckT_relErr_vec_mat.txt")


##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.plot(Ttrue, "o", linewidth = 2, label = "Ttrue")
plt.plot(recT, "o", linewidth = 2, label = "recT")
plt.plot(toCheckT, "o", linewidth = 2, label = "toCheckT")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Temperature [K]', fontsize=25)
#plt.xlim(0,5)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig = plt.figure(2,figsize=(12,8))
plt.semilogy(recT_relErr, "o", linewidth = 2, label = "recT")
plt.semilogy(toCheckT_relErr, "o", linewidth = 2, label = "toCheckT")
plt.semilogy(toCheckTrecT_relErr, "o", linewidth = 2, label = "toCheckTrecT")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Relative Errors []', fontsize=25)
#plt.xlim(0,5)
plt.grid()
plt.legend()


##############################################################
plt.show()

