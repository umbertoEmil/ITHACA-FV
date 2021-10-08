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
probeCenter_gRec = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_gRec_mat.txt")
probeCenter_gTrue = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_gTrue_mat.txt")
probeCenter_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_Ttrue_mat.txt")
probeCenter_Trec = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_Trec_mat.txt")
probe2_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/probe2_Ttrue_mat.txt")
probe2_Trec = np.loadtxt("./ITHACAoutput/testInverse/probe2_Trec_mat.txt")

probeCenter_gRec_noise01 = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_gRec_noise01_mat.txt")
probeCenter_Trec_noise01 = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_Trec_noise01_mat.txt")
probe2_Trec_noise01 = np.loadtxt("./ITHACAoutput/testInverse/probe2_Trec_noise01_mat.txt")


probeCenter_gRec_noise001 = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_gRec_noise001_mat.txt")
probeCenter_Trec_noise001 = np.loadtxt("./ITHACAoutput/testInverse/probeCenter_Trec_noise001_mat.txt")
probe2_Trec_noise001 = np.loadtxt("./ITHACAoutput/testInverse/probe2_Trec_noise001_mat.txt")

#relErr_L2 = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")
#relErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/relErrLinfNorm_mat.txt")
#avgRelErr_L2 =   np.loadtxt("./ITHACAoutput/testInverse/avgRelErrL2norm_mat.txt")
#avgRelErr_Linf = np.loadtxt("./ITHACAoutput/testInverse/avgRelErrLinfNorm_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.plot(t, probeCenter_gTrue, "b", linewidth = 2, label = r"$g_t$")
#plt.plot(t, probeCenter_gRec, "k--", linewidth = 2, label = r"$g$")
plt.plot(t, probeCenter_gRec_noise001, "C2--",  linewidth = 2, label = r"$g(\zeta = 1e-4)$")
plt.plot(t, probeCenter_gRec_noise01, "C1-.", linewidth = 2, label = r"$g(\zeta = 1e-3)$")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel(r'Heat flux [$W/m^2$]', fontsize=25)
plt.xlim(0,50)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

##############################################################
fig, axes = plt.subplots(figsize=(12,8))

plt.plot(t, probeCenter_Ttrue, "bo-", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None",  label = "at (1, 0, 0.6)")
#plt.plot(t, probeCenter_Trec, "bs--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, probeCenter_Trec_noise001, "C2d--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, probeCenter_Trec_noise01, "C1X-.", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
plt.plot(t, probe2_Ttrue, "ko-", markevery=10, linewidth = 2, label = "at (1, 0.02, 0.6)", markersize=8, markerfacecolor="None")
#plt.plot(t, probe2_Trec, "ks--", markevery=10, linewidth = 2, markersize=8, markerfacecolor="None")
plt.plot(t, probe2_Trec_noise001,"C2d--", markevery=10, linewidth = 2, markersize=8, markerfacecolor="None")
plt.plot(t, probe2_Trec_noise01, "C1X-.",markevery=10, linewidth = 2, markersize=8, markerfacecolor="None")


plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Temperature [K]', fontsize=25)
plt.xlim(0,50)
plt.grid()

leg = plt.legend()
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle=i, marker=j, linewidth = 2, markerfacecolor="None", ls="")[0] for i in ["-", "--", "-."] for j in ["o", "d", "X"]]
plt.legend(handles=h, labels=["True", r"$\zeta = 1e-4$", r"$\zeta = 1e-3$"], ncol = 3, bbox_to_anchor=(0.05, 0.08),loc=2, borderaxespad=0.)

##############################################################
##fig = plt.figure(3,figsize=(12,8))
##plt.semilogy(t, relErr_L2, "kx--", linewidth = 2, markevery=10, label = "FOM")
##plt.semilogy(t, relErr_Linf, "ko", linewidth = 2, markevery=10)
##
##leg = plt.legend()
###axes.add_artist(leg)
###h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
###plt.legend(handles=h, labels=["Linf", "L2"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
##plt.xlabel('Time [s]', fontsize=25)
##plt.grid()
##
################################################################
##fig = plt.figure(5,figsize=(12,8))
##plt.semilogy(t, avgRelErr_L2, "k--", linewidth = 2, label = "FOM")
##plt.semilogy(t, avgRelErr_Linf, "k", linewidth = 2)
##
##leg = plt.legend()
###axes.add_artist(leg)
###h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
###plt.legend(handles=h, labels=["Linf", "L2"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
##plt.xlabel('Time [s]', fontsize=25)
##plt.grid()

##############################################################
plt.show()

