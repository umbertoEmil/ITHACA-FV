import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np
import sys
sys.path.insert(0, "./")

import scipy.stats

import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc
 
def plot_mean_and_CI(mean, lb, ub, color_mean=None, color_shading=None):
    # plot the shaded range of the confidence intervals
    plt.fill_between(range(mean.shape[0]), ub, lb,
                     color=color_shading, alpha=.5)
    # plot the mean on top
    plt.plot(mean, color_mean)



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
forcer = np.loadtxt("./ITHACAoutput/forcer_mat.txt")
stateRec = np.loadtxt("./ITHACAoutput/stateRec_mat.txt")
singleVariableSample = np.loadtxt("./ITHACAoutput/singleVariableSamples_mat.txt")
posteriorMean = np.loadtxt("./ITHACAoutput/posteriorMean_mat.txt")
minConfidence = np.loadtxt("./ITHACAoutput/minConfidence_mat.txt")
maxConfidence = np.loadtxt("./ITHACAoutput/maxConfidence_mat.txt")


fig = plt.figure(1,figsize=(8,6))
plt.plot(time, X[0,:], linewidth = 2, label="1")
plt.plot(time, X[1,:], linewidth = 2, label="2")
plt.plot(time, X[2,:], linewidth = 2, label="3")
plt.plot(time, X[3,:], linewidth = 2, label="4")

plt.plot(time, stateRec[0,:], "o", linewidth = 2, label="p1")
plt.plot(time, stateRec[1,:], "o", linewidth = 2, label="p2")
plt.plot(time, stateRec[2,:], "o", linewidth = 2, label="p3")
plt.plot(time, stateRec[3,:], "o", linewidth = 2, label="p4")

plt.legend()
plt.xlabel('Time [s]', fontsize=25)
plt.grid()

fig = plt.figure(2,figsize=(8,6))
plt.grid()

plt.plot(time, forcer[0,:], linewidth = 2, label="1")
plt.plot(time, forcer[1,:], linewidth = 2, label="2")
plt.plot(time, forcer[2,:], linewidth = 2, label="3")
plt.plot(time, forcer[3,:], linewidth = 2, label="4")
plt.plot(time, posteriorMean[0,:], "o", linewidth = 2, label="p1")
plt.plot(time, posteriorMean[1,:], "o", linewidth = 2, label="p2")
plt.plot(time, posteriorMean[2,:], "o", linewidth = 2, label="p3")
plt.plot(time, posteriorMean[3,:], "o", linewidth = 2, label="p4")
plt.legend()

fig = plt.figure(3,figsize=(8,6))
for column in singleVariableSample.T:
    plt.plot(time, column, "ko", markersize=0.8)
plt.fill_between(time, minConfidence[0,:], maxConfidence[0,:], color='b', alpha=.1, label= "b(1)")
plt.plot(time,posteriorMean[0,:], color='b')

plt.legend()

plt.show()





















