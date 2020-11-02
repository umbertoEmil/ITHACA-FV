import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)




cellSize = np.array([0.25, 0.1, 0.0667, 0.05, 0.03333, 0.025])



################################################
############## inverseProblem ##################
################################################
inverseError_L2norm = np.array([
    0.0067407202887057532728,
    0.042714151965948153611,
    0.033608975582644841362,
    0.032302701782013193421,
    0.028675342920687824089,
    0.028570443774623317801
])



inverseError_LinfNorm = np.array([
    0.0089291945991111440484,
    0.13466349985925932242,
    0.08271849056480282125,
    0.087755049438383886384,
    0.075432247804407520642,
    0.080984854045179119342
])

#inverseError_L2norm = np.array([
#0.03892633138862037,
#0.03892633138862037,
#0.027343409995907095,
#0.02502546453939964,
#0.02379717836951923,
#0.023323013002922565
#])
#
#inverseError_LinfNorm = np.array([
#0.1677838040888888,
#0.1677838040888888,
#0.09589284076625712,
#0.06873519795604392,
#0.06421078016382116,
#0.06568411651570251
#])


f = plt.figure(20,figsize=(12,8))
ax1 = plt.axes(xscale='log')
ax1.plot(cellSize, inverseError_L2norm, 'bo', markersize=15, label=r'$||\epsilon||_{L^2(\Gamma_{in})}$')
ax1.plot(cellSize, inverseError_LinfNorm, 'kv', markersize=15, label=r'$||\epsilon||_{L^\infty(\Gamma_{in})}$')
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
plt.xticks(cellSize)
plt.xlabel("Cell edge lenght [m]", fontsize=25)
plt.legend(loc='best', fontsize=25)
plt.grid()




plt.show()