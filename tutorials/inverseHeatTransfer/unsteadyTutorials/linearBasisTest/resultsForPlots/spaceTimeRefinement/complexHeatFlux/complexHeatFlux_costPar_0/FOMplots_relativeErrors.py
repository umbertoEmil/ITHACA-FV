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

deltaX = "30 10 25"
deltaT = ["0.1", "0.2", "0.25", "0.5"]

i = 1

t_final = 50

for dt in deltaT:
    dT = float(dt)

    heatFluxRelErr_L2 = np.loadtxt("./heatFluxRelErr_L2_" + deltaX + "_" + dt)
    t = np.linspace(dT, t_final, heatFluxRelErr_L2.size, endpoint=True)
    print t
    probeCenter_trueHeatFlux = np.loadtxt("./probeHeat_true_" + dt)
    probeCenter_gRec = np.loadtxt("./probeHeat_rec_" + deltaX + "_" + dt)

    fig = plt.figure(i,figsize=(12,8))
    plt.plot(t, probeCenter_trueHeatFlux, "r", linewidth = 2, label = r"$g_t$")
    plt.plot(t, probeCenter_gRec, "bo-", linewidth = 2, label = r"$g_{FOM}$")
    
    
    plt.xlabel('Time [s]', fontsize=25)
    plt.ylabel(r'Heat flux $[W/m^2]$', fontsize=25)
    plt.xlim(0,50)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.grid()
    plt.legend()

    fig = plt.figure(i + 100,figsize=(12,8))
    plt.semilogy(t,heatFluxRelErr_L2,   "bo-", linewidth = 2,markersize = 6)
    
    
    plt.xlabel('Time [s]', fontsize=25)
    plt.ylabel(r"$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$", fontsize=25)
    plt.xlim(0,50)
    plt.grid()
    i = i + 1

#
#t = np.loadtxt("./ITHACAoutput/testInverse/FOM/timeSteps_mat.txt")
#probeCenter_gRec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_gRec_mat.txt")
#probeCenter_trueHeatFlux = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_trueHeatFlux_mat.txt")
#probeCenter_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_Ttrue_mat.txt")
#probeCenter_Trec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probeCenter_Trec_mat.txt")
#probe2_Ttrue = np.loadtxt("./ITHACAoutput/testInverse/FOM/probe2_Ttrue_mat.txt")
#probe2_Trec = np.loadtxt("./ITHACAoutput/testInverse/FOM/probe2_Trec_mat.txt")
#onlineWindowsVec = np.loadtxt("./ITHACAoutput/testInverse/ROM/onlineWindowsVec_mat.txt")
#offlineWindowsVec = np.loadtxt("./ITHACAoutput/testInverse/ROM/offlineWindowsVec_mat.txt")
#
#REDprobeCenter_gRec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probeCenter_gRec_mat.txt")
#REDprobeCenter_Trec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probeCenter_Trec_mat.txt")
#REDprobe2_Trec = np.loadtxt( "./ITHACAoutput/testInverse/ROM/probe2_Trec_mat.txt")
#
#heatFluxRelErr_L2_FOM = np.loadtxt("./ITHACAoutput/testInverse/FOM/heatFluxRelErr_L2_mat.txt")
#heatFluxRelErr_Linf_FOM = np.loadtxt("./ITHACAoutput/testInverse/FOM/heatFluxRelErr_Linf_mat.txt")
#
#heatFluxRelErr_L2_ROM = np.loadtxt("./ITHACAoutput/testInverse/ROM/heatFluxRelErr_L2_mat.txt")
#heatFluxRelErr_Linf_ROM = np.loadtxt("./ITHACAoutput/testInverse/ROM/heatFluxRelErr_Linf_mat.txt")
#
#timeVec = np.loadtxt("./ITHACAoutput/testInverse/ROM/timeVec_mat.txt")
#
###############################################################
#fig = plt.figure(2,figsize=(12,8))
#plt.plot(t, probeCenter_trueHeatFlux, "r", linewidth = 2, label = r"$g_t$")
#plt.plot(t, probeCenter_gRec, "bo-", linewidth = 2, label = r"$g_{FOM}$")
#plt.plot(t, REDprobeCenter_gRec, "ks--", linewidth = 2, label = r"$g_{ROM}$")
#
#
#plt.xlabel('Time [s]', fontsize=25)
#plt.ylabel(r'Heat flux $[W/m^2]$', fontsize=25)
#plt.xlim(0,50)
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.grid()
#plt.legend()
#
#deltaT = 1 #size of time windows
#for x in offlineWindowsVec:
#    plt.axvline(x = x)
#    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.1)
#
###############################################################
#fig, axes =  plt.subplots(figsize=(12,8))
#plt.semilogy(timeVec,heatFluxRelErr_L2_FOM,   "bo-", linewidth = 2,markersize = 6)
#plt.semilogy(timeVec,heatFluxRelErr_L2_FOM,   "b-", linewidth = 2,markersize = 6, label = "FOM")
#plt.semilogy(timeVec,heatFluxRelErr_Linf_FOM, "bs-", linewidth = 2,markersize = 6)
#plt.semilogy(timeVec,heatFluxRelErr_L2_ROM,   "ko--",   linewidth = 2,markersize = 6)
#plt.semilogy(timeVec,heatFluxRelErr_L2_ROM,   "k--",   linewidth = 2,markersize = 6, label = "ROM")
#plt.semilogy(timeVec,heatFluxRelErr_Linf_ROM, "ks--", linewidth = 2,markersize = 6)
#
#plt.xlabel('Time [s]', fontsize=25)
#plt.xlim(0,50)
##plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.grid()
#
#leg = plt.legend(loc='best')
#axes.add_artist(leg)
#h = [plt.plot([],[], color="k", linestyle="-", marker=j, linewidth = 2, markerfacecolor="k", markersize = 15, ls="")[0] for j in ["o", "s"]]
##plt.legend(handles=h, labels=["Piecewise constant", "Piecewise linear"], ncol = 2, bbox_to_anchor=(0.001, 0.06),loc=2, borderaxespad=0.)
#plt.legend(handles=h, labels=[r"$||e_{rel}||_{L^2(\Gamma_{s_{in}})}$", r"$||e_{rel}||_{L^\infty(\Gamma_{s_{in}})}$"], ncol = 2, bbox_to_anchor=(0., 1.056), frameon=False, edgecolor = ('k'), loc="upper left" , borderaxespad=0.)
#
#deltaT = 1 #size of time windows
#for x in offlineWindowsVec:
#    plt.axvline(x = x)
#    plt.axvspan(x - deltaT, x, facecolor='#DC143C', alpha=0.1)

###############################################################
#fig, axes = plt.subplots(figsize=(12,8))
#
#plt.plot(t, probeCenter_Ttrue, "bo-", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None", label = "at (1, 0, 0.6)")
#plt.plot(t, probeCenter_Trec, "bs--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
#plt.plot(t, REDprobeCenter_Trec, "rX-.", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
#plt.plot(t, probe2_Ttrue, "ko-", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None",  label = "at (1, 0.015, 0.6)")
#plt.plot(t, probe2_Trec, "ks--", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
#plt.plot(t, REDprobe2_Trec, "rX-.", linewidth = 2, markevery=10, markersize=8, markerfacecolor="None")
#
#plt.xlabel('Time [s]', fontsize=25)
#plt.xlim(0,50)
#plt.ylabel('Temperature [K]', fontsize=25)
#plt.grid()
#
#
#leg = plt.legend()
#axes.add_artist(leg)
##h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
#h = [plt.plot([],[], color="k", linestyle=i, marker=j, linewidth = 2, markerfacecolor="None", ls="")[0] for i in ["-", "--", "-."] for j in ["o", "s", "X"]]
#plt.legend(handles=h, labels=["True", "Full", "Reduced"], ncol = 3, bbox_to_anchor=(0.05, 0.08),loc=2, borderaxespad=0.)

plt.show()

