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

time = np.loadtxt("./ITHACAoutput/time_mat.txt")
X = np.loadtxt("./ITHACAoutput/X_mat.txt")
posteriorMean = np.loadtxt("./ITHACAoutput/posteriorMean_mat.txt")


fig = plt.figure(1,figsize=(8,6))
plt.plot(time, X[0,:], linewidth = 2, label="1")
plt.plot(time, X[1,:], linewidth = 2, label="2")
plt.plot(time, X[2,:], linewidth = 2, label="3")
plt.plot(time, X[3,:], linewidth = 2, label="4")

plt.plot(time, posteriorMean[0,:], "o", linewidth = 2, label="p1")
plt.plot(time, posteriorMean[1,:], "o", linewidth = 2, label="p2")
plt.plot(time, posteriorMean[2,:], "o", linewidth = 2, label="p3")
plt.plot(time, posteriorMean[3,:], "o", linewidth = 2, label="p4")

plt.legend()
plt.xlabel('Time [s]', fontsize=25)
plt.grid()

plt.show()
