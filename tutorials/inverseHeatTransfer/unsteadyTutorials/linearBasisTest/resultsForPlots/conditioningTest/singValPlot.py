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
delta_x = ["100 20 85", "50 20 45", "40 15 35", "30 10 25", "20 5 15"]
deltaT = ["0.1", "0.2", "0.25", "0.5"]

fig = plt.figure(1,figsize=(12,8))
for dx in delta_x:
    condNum = []
    for dt in deltaT:
        condNum.append(np.loadtxt("./conditioningNumber_" + dx + "_" + dt))
    plt.semilogy(deltaT, condNum, "o-",linewidth = 2, markersize=15, label = r"$Mesh = $" + dx)

plt.xlabel(r"$\Delta t$", fontsize=25)
plt.ylabel(r"Conditioning number = $\frac{\sigma_{max}}{\sigma_{min}}$", fontsize=25)
plt.grid()
plt.legend()



fig = plt.figure(2,figsize=(12,8))
for dt in deltaT:
    singVal = np.loadtxt("./singularValues_" + deltaX + "_" + dt)
    plt.semilogy(singVal, "o",linewidth = 2, label = r"$\Delta t = $" + dt)
    
    plt.ylabel("Singular Values", fontsize=25)
plt.grid()
plt.legend()

##############################################################
plt.show()

