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
minConfidence = np.loadtxt("./ITHACAoutput/minConfidence_mat.txt")
maxConfidence = np.loadtxt("./ITHACAoutput/maxConfidence_mat.txt")



fig = plt.figure(1,figsize=(8,6))
plt.plot(time, X[0,:],"b--", linewidth = 2, label="1")
plt.plot(time, X[1,:],"k--", linewidth = 2, label="2")
plt.plot(time, X[2,:],"g--", linewidth = 2, label="3")
plt.plot(time, X[3,:],"m--", linewidth = 2, label="4")

#plt.plot(time, posteriorMean[0,:], "o", linewidth = 2, label="p1")
#plt.plot(time, posteriorMean[1,:], "o", linewidth = 2, label="p2")
#plt.plot(time, posteriorMean[2,:], "o", linewidth = 2, label="p3")
#plt.plot(time, posteriorMean[3,:], "o", linewidth = 2, label="p4")

plt.fill_between(time, minConfidence[0,:], maxConfidence[0,:], color='b', alpha=.1, label="xrec1" )
plt.plot(time,posteriorMean[0,:], color='b')
plt.fill_between(time, minConfidence[1,:], maxConfidence[1,:], color='k', alpha=.1, label="xrec2" )
plt.plot(time,posteriorMean[1,:], color='k')
plt.fill_between(time, minConfidence[2,:], maxConfidence[2,:], color='g', alpha=.1, label="xrec3" )
plt.plot(time,posteriorMean[2,:], color='g')
plt.fill_between(time, minConfidence[2,:], maxConfidence[2,:], color='m', alpha=.1, label="xrec3" )
plt.plot(time,posteriorMean[3,:], color='m')

plt.legend()
plt.xlabel('Time [s]', fontsize=25)
plt.grid()


plt.show()
