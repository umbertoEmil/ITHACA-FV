import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
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

time_C = np.loadtxt("./coarseMesh/ITHACAoutput/test/timeVec_mat.txt")
relErr_L2 = np.loadtxt("./coarseMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf = np.loadtxt("./coarseMesh/ITHACAoutput/test/relErr_Linf_mat.txt")

time_M = np.loadtxt("./mediumMesh/ITHACAoutput/test/timeVec_mat.txt")
relErr_L2_M = np.loadtxt("./mediumMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf_M = np.loadtxt("./mediumMesh/ITHACAoutput/test/relErr_Linf_mat.txt")

time_F = np.loadtxt("./fineMesh/ITHACAoutput/test/timeVec_mat.txt")
relErr_L2_F = np.loadtxt("./fineMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf_F = np.loadtxt("./fineMesh/ITHACAoutput/test/relErr_Linf_mat.txt")




fig, axes = plt.subplots()
plt.semilogy(time_C, relErr_L2,"b", linewidth = 2)
plt.semilogy(time_C, relErr_Linf,"k", linewidth = 2, label=r"$\Delta t = 1s$")
plt.semilogy(time_M, relErr_L2_M,"b--", linewidth = 2)
plt.semilogy(time_M, relErr_Linf_M,"k--", linewidth = 2, label=r"$\Delta t = 0.1s$")
plt.semilogy(time_F, relErr_L2_F,"b.", linewidth = 2)
plt.semilogy(time_F, relErr_Linf_F,"k.", linewidth = 2, label=r"$\Delta t = 0.01s$")

leg = plt.legend(loc="upper left", fontsize=25)
axes.add_artist(leg)
h = [plt.plot([],[], color=i, marker='o', markersize=15, ls="")[0] for i in ["b", "k"]]# for j in ["-" "--"]]
plt.legend(fontsize=25, handles=h, labels=[r"$||\epsilon||_{L^2(\Omega_s)}$", r"$||\epsilon||_{L^\infty(\Omega_s)}$"], ncol = 3, bbox_to_anchor=(0., 1.11),loc=2, borderaxespad=0.)


#axes.set_title('Implicit Euler', fontsize=25, y=.05, x = 0.7)
plt.xlabel('Time [s]', fontsize=25)
plt.xlim(0,time_C[-1])
#plt.title("CrankNicolson 0.9", fontsize=25)
plt.grid()

plt.show()
