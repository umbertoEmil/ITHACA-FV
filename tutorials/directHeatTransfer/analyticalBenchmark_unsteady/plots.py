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

time = np.loadtxt("./coarseMesh/ITHACAoutput/test/timeVec_mat.txt")
relErr_L2 = np.loadtxt("./coarseMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf = np.loadtxt("./coarseMesh/ITHACAoutput/test/relErr_Linf_mat.txt")
relErr_L2_M = np.loadtxt("./mediumMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf_M = np.loadtxt("./mediumMesh/ITHACAoutput/test/relErr_Linf_mat.txt")
relErr_L2_F = np.loadtxt("./fineMesh/ITHACAoutput/test/relErr_L2_mat.txt")
relErr_Linf_F = np.loadtxt("./fineMesh/ITHACAoutput/test/relErr_Linf_mat.txt")




fig = plt.figure(1,figsize=(8,6))
plt.semilogy(time, relErr_L2,"b", linewidth = 2, label="L2")
plt.semilogy(time, relErr_Linf,"k", linewidth = 2, label="Linf")
plt.semilogy(time, relErr_L2_M,"b--", linewidth = 2)
plt.semilogy(time, relErr_Linf_M,"k--", linewidth = 2)
plt.semilogy(time, relErr_L2_F,"b.", linewidth = 2)
plt.semilogy(time, relErr_Linf_F,"k.", linewidth = 2)

plt.xlabel('Time [s]', fontsize=25)
plt.grid()
plt.legend()

plt.show()
