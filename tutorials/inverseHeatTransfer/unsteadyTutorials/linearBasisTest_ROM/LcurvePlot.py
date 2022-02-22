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

eta30 = np.loadtxt("./ITHACAoutput/Lcurve/eta30_mat.txt")
rho30 = np.loadtxt("./ITHACAoutput/Lcurve/rho30_mat.txt")
eta20 = np.loadtxt("./ITHACAoutput/Lcurve/eta20_mat.txt")
rho20 = np.loadtxt("./ITHACAoutput/Lcurve/rho20_mat.txt")
eta10 = np.loadtxt("./ITHACAoutput/Lcurve/eta10_mat.txt")
rho10 = np.loadtxt("./ITHACAoutput/Lcurve/rho10_mat.txt")
regParam = np.loadtxt("./ITHACAoutput/Lcurve/regParam_mat.txt")

##############################################################
fig = plt.figure(1,figsize=(12,8))
plt.loglog(rho30,eta30, linewidth = 2, label = "30")
plt.loglog(rho20,eta20, linewidth = 2, label = "20")
plt.loglog(rho10,eta10, linewidth = 2, label = "10")

plt.xlabel('rho', fontsize=25)
plt.ylabel('eta', fontsize=25)
plt.grid()
plt.legend()



##############################################################
plt.show()

