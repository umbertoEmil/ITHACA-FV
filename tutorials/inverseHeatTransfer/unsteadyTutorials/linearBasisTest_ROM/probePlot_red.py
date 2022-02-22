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

t = np.loadtxt("./ITHACAoutput/testInverse/FOM/timeSteps_mat.txt")
probeCenter_gRec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_gRec_mat.txt")
probeCenter_trueHeatFlux = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_trueHeatFlux_mat.txt")
probeCenter_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_Ttrue_mat.txt")
probeCenter_Trec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_Trec_mat.txt")
probe2_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/FOM/probe2_Ttrue_mat.txt")
probe2_Trec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probe2_Trec_mat.txt")
onlineWindowsVec = np.loadtxt("./ITHACAoutput/testInverse/ROM/onlineWindowsVec_mat.txt")
offlineWindowsVec = np.loadtxt("./ITHACAoutput/testInverse/ROM/offlineWindowsVec_mat.txt")

REDprobeCenter_gRec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probeCenter_gRec_mat.txt")
REDprobeCenter_Trec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probeCenter_Trec_mat.txt")
REDprobe2_Trec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probe2_Trec_mat.txt")


#relErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")
#relErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/relErrLinfNorm_mat.txt")
#avgRelErr_L2 =   np.loadtxt("./ITHACAoutput/testInverse/avgRelErrL2norm_mat.txt")
#avgRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/avgRelErrLinfNorm_mat.txt")
#REDrelErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/ROM/relErrL2norm_mat.txt")
#REDrelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/ROM/relErrLinfNorm_mat.txt")
#REDavgRelErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/ROM/avgRelErrL2norm_mat.txt")
#REDavgRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/ROM/avgRelErrLinfNorm_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.plot(t, probeCenter_trueHeatFlux, "b", linewidth = 2, label = r"$g_t$")
plt.plot(t, probeCenter_gRec, "k--", linewidth = 2, label = "Full")
plt.plot(t, REDprobeCenter_gRec, "r--", linewidth = 2, label = "Reduced")


plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux $[W/m^2]$', fontsize=25)
plt.xlim(0,50)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

deltaT = 1 #size of time windows
for x in offlineWindowsVec:
    plt.axvline(x = x)
    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.1)

##############################################################
fig, axes = plt.subplots(figsize=(12,8))

plt.plot(t, probeCenter_Ttrue, "bo-", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None", label = "at (1, 0, 0.6)")
plt.plot(t, probeCenter_Trec, "bs--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, REDprobeCenter_Trec, "rX-.", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, probe2_Ttrue, "ko-", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None",  label = "at (1, 0.015, 0.6)")
plt.plot(t, probe2_Trec, "ks--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, REDprobe2_Trec, "rX-.", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")

plt.xlabel('Time [s]', fontsize=25)
plt.xlim(0,50)
plt.ylabel('Temperature [K]', fontsize=25)
plt.grid()


leg = plt.legend()
axes.add_artist(leg)
#h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
h = [plt.plot([],[], color="k", linestyle=i, marker=j, linewidth = 2, markerfacecolor="None", ls="")[0] for i in ["-", "--", "-."] for j in ["o", "s", "X"]]
plt.legend(handles=h, labels=["True", "Full", "Reduced"], ncol = 3, bbox_to_anchor=(0.05, 0.08),loc=2, borderaxespad=0.)


###############################################################
#fig = plt.figure(3,figsize=(12,8))
#plt.semilogy(t, relErr_L2, "k--", linewidth = 2, label = "FOM")
#plt.semilogy(t, REDrelErr_L2, "r--", linewidth = 2, label = "ROM")
#plt.semilogy(t, relErr_Linf, "k", linewidth = 2)
#plt.semilogy(t, REDrelErr_Linf, "r", linewidth = 2)
#
#leg = plt.legend()
##axes.add_artist(leg)
##h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
##plt.legend(handles=h, labels=["Linf", "L2"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
#plt.xlabel('Time [s]', fontsize=25)
#plt.grid()
#
#for x in offlineWindowsVec:
#    plt.axvline(x = x)
#    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.2)
#
###############################################################
#fig = plt.figure(4,figsize=(12,8))
#REDErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/ROM/REDErrL2norm_mat.txt")
#REDErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/ROM/REDErrLinfNorm_mat.txt")
#
#plt.semilogy(t, REDErr_L2, "k--", linewidth = 2, label = "L2")
#plt.semilogy(t, REDErr_Linf, "k", linewidth = 2, label = "Linf")
#
#leg = plt.legend()
#plt.xlabel('Time [s]', fontsize=25)
#plt.grid()
#
#
#for x in offlineWindowsVec:
#    plt.axvline(x = x)
#    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.2)
#
###############################################################
#fig = plt.figure(5,figsize=(12,8))
#plt.semilogy(t, avgRelErr_L2, "k--", linewidth = 2, label = "FOM")
#plt.semilogy(t, REDavgRelErr_L2, "r--", linewidth = 2, label = "ROM")
#plt.semilogy(t, avgRelErr_Linf, "k", linewidth = 2)
#plt.semilogy(t, REDavgRelErr_Linf, "r", linewidth = 2)
#
#leg = plt.legend()
##axes.add_artist(leg)
##h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
##plt.legend(handles=h, labels=["Linf", "L2"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
#plt.xlabel('Time [s]', fontsize=25)
#plt.grid()
#
#for x in offlineWindowsVec:
#    plt.axvline(x = x)
#    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.2)

##############################################################
plt.show()

