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
probe1_gRec = np.loadtxt("./ITHACAoutput/testInverse/probe1_gRec_mat.txt")
probe1_gTrue = np.loadtxt("./ITHACAoutput/testInverse/probe1_gTrue_mat.txt")
probe1_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/probe1_Ttrue_mat.txt")
probe1_Trec = np.loadtxt("./ITHACAoutput/testInverse/probe1_Trec_mat.txt")
probe2_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/probe2_Ttrue_mat.txt")
probe2_Trec = np.loadtxt("./ITHACAoutput/testInverse/probe2_Trec_mat.txt")
onlineWindowsVec = np.loadtxt("./ITHACAoutput/testInverseReduced/onlineWindowsVec_mat.txt")

REDprobe1_gRec = np.loadtxt( "./ITHACAoutput/testInverseReduced/probe1_gRec_mat.txt")
REDprobe1_Trec = np.loadtxt( "./ITHACAoutput/testInverseReduced/probe1_Trec_mat.txt")
REDprobe2_Trec = np.loadtxt( "./ITHACAoutput/testInverseReduced/probe2_Trec_mat.txt")


relErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")
relErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/relErrLinfNorm_mat.txt")
REDrelErr_L2 = np.loadtxt("./ITHACAoutput/testInverseReduced/relErrL2norm_mat.txt")
REDrelErr_Linf = np.loadtxt("./ITHACAoutput/testInverseReduced/relErrLinfNorm_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.plot(t, probe1_gTrue, "b", linewidth = 2, label = "True")
plt.plot(t, probe1_gRec, "k--", linewidth = 2, label = "Estimated")
plt.plot(t, REDprobe1_gRec, "r--", linewidth = 2, label = "RED")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Heat flux [W/m2]', fontsize=25)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

deltaT = 1 #size of time windows
for x in onlineWindowsVec:
    plt.axvline(x = x)
    plt.axvspan(x - deltaT, x, facecolor='#2ca02c', alpha=0.2)

##############################################################
fig, axes = plt.subplots(figsize=(12,8))

plt.plot(t, probe1_Ttrue, "b", linewidth = 2, label = "at (0, 0, 0.6)")
plt.plot(t, probe1_Trec, "b--", linewidth = 2)
plt.plot(t, REDprobe1_Trec, "r--", linewidth = 2)
plt.plot(t, probe2_Ttrue, "k", linewidth = 2, label = "at (0, 0.015, 0.6)")
plt.plot(t, probe2_Trec, "k--", linewidth = 2)
plt.plot(t, REDprobe2_Trec, "c--", linewidth = 2)

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Temperature [K]', fontsize=25)
plt.grid()

leg = plt.legend()
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
plt.legend(handles=h, labels=["True", "Estimated"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
for x in onlineWindowsVec:
    plt.axvline(x = x)
    plt.axvspan(x - deltaT, x, facecolor='#2ca02c', alpha=0.2)

##############################################################
fig = plt.figure(3,figsize=(12,8))
plt.semilogy(t, relErr_L2, "k--", linewidth = 2, label = "FOM")
plt.semilogy(t, REDrelErr_L2, "r--", linewidth = 2, label = "ROM")
plt.semilogy(t, relErr_Linf, "k", linewidth = 2)
plt.semilogy(t, REDrelErr_Linf, "r", linewidth = 2)

leg = plt.legend()
#axes.add_artist(leg)
#h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
#plt.legend(handles=h, labels=["Linf", "L2"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
plt.xlabel('Time [s]', fontsize=25)
plt.grid()
for x in onlineWindowsVec:
    plt.axvline(x = x)
    plt.axvspan(x - deltaT, x, facecolor='#2ca02c', alpha=0.2)

##############################################################
fig = plt.figure(4,figsize=(12,8))
REDErr_L2 = np.loadtxt("./ITHACAoutput/testInverseReduced/REDErrL2norm_mat.txt")
REDErr_Linf = np.loadtxt("./ITHACAoutput/testInverseReduced/REDErrLinfNorm_mat.txt")

plt.semilogy(t, REDErr_L2, "k--", linewidth = 2, label = "L2")
plt.semilogy(t, REDErr_Linf, "k", linewidth = 2, label = "Linf")

leg = plt.legend()
plt.xlabel('Time [s]', fontsize=25)
plt.grid()


##############################################################
plt.show()

