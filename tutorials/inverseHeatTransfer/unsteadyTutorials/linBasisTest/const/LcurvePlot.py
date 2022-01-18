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

eta0 = np.loadtxt("./ITHACAoutput/Lcurve/eta0_mat.txt")
rho0 = np.loadtxt("./ITHACAoutput/Lcurve/rho0_mat.txt")
eta2 = np.loadtxt("./ITHACAoutput/Lcurve/eta2_mat.txt")
rho2 = np.loadtxt("./ITHACAoutput/Lcurve/rho2_mat.txt")
eta1 = np.loadtxt("./ITHACAoutput/Lcurve/eta1_mat.txt")
rho1 = np.loadtxt("./ITHACAoutput/Lcurve/rho1_mat.txt")
regParam = np.loadtxt("./ITHACAoutput/Lcurve/regParam_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.loglog(rho0,eta0, linewidth = 2, label = "0")
plt.loglog(rho1,eta1, linewidth = 2, label = "1")
plt.loglog(rho2,eta2, linewidth = 2, label = "2")

plt.xlabel('rho', fontsize=25)
plt.ylabel('eta', fontsize=25)
plt.grid()
plt.legend()



##############################################################
plt.show()

