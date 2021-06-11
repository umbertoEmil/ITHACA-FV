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



AbsError_L2norm = np.array([
6.18122e-07,
3.64878e-07,
3.19016e-07,
2.93347e-07
])

AbsError_LinfNorm = np.array([
2.17167e-06,
1.28802e-06,
1.15261e-06,
1.06198e-06
])

gridSize = np.array([20**3, 30**3, 40**3, 50**3])

f = plt.figure(1,figsize=(12,8))
ax1 = plt.axes(xscale='log')
ax1.plot(gridSize,AbsError_L2norm, 'ko', markersize=15, label=r'$L^2$ norm')
ax1.plot(gridSize,AbsError_LinfNorm, 'bv', markersize=15, label=r'$L^\infty$ norm')
ax1.set_xticks(gridSize)
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
ax1.set_xlabel("Mesh size", fontsize=25)
leg = plt.legend(loc='best', fontsize=25)


#plt.semilogx(gridSize,AbsError_L2norm, 'bo')
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.xticks(gridSize)
plt.grid()
plt.show()
