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

## STEADY
probe1_gRec_steady = np.loadtxt("./ITHACAoutput/steadyTest/probe1_gRec_mat.txt")
probe1_Trec_steady = np.loadtxt("./ITHACAoutput/steadyTest/probe1_Trec_mat.txt")
probe2_Trec_steady = np.loadtxt("./ITHACAoutput/steadyTest/probe2_Trec_mat.txt")
time_steady = np.loadtxt("./ITHACAoutput/steadyTest/time_mat.txt")

fig = plt.figure(1,figsize=(12,8))
plt.plot(t, probe1_gTrue, "b", linewidth = 2, label = "True")
plt.plot(t, probe1_gRec, "k--", linewidth = 2, label = "Estimated")
plt.plot(time_steady, probe1_gRec_steady, "mo", label = "Steady")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Heat flux [W/m2]', fontsize=25)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()

#fig = plt.figure(10,figsize=(8,6))
fig, axes = plt.subplots(figsize=(12,8))

plt.plot(t, probe1_Ttrue, "b", linewidth = 2, label = "at (0, 0, 0.6)")
plt.plot(t, probe1_Trec, "b--", linewidth = 2)
plt.plot(t, probe2_Ttrue, "k", linewidth = 2, label = "at (0, 0.015, 0.6)")
plt.plot(t, probe2_Trec, "k--", linewidth = 2)
plt.plot(time_steady, probe1_Trec_steady, "mo", label = "Steady")
plt.plot(time_steady, probe2_Trec_steady, "mo")

plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('Temperature [K]', fontsize=25)
plt.grid()


leg = plt.legend()
axes.add_artist(leg)
h = [plt.plot([],[], color="k", linestyle=i, linewidth = 2, ls="")[0] for i in ["-", "--"]]# for j in ["-" "--"]]
plt.legend(handles=h, labels=["True", "Estimated"], ncol = 2, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)

plt.show()

