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

matrix = np.loadtxt('ITHACAoutput/MHMCMCregParTestResults/Js_mat.txt')

fig, ax = plt.subplots(figsize=(10, 6))
fig.subplots_adjust(bottom=0.2, left=0.2)
#ax1 = plt.axes(xscale='log')
#ax1.plot(matrix[:,0], matrix[:,1], 'k', linewidth=2)
#ax1.plot(gridSize,AbsError_LinfNorm, 'bv', markersize=15, label=r'$L^\infty$ norm')
#ax1.set_xticks(gridSize)
#ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#ax1.set_xlabel("Mesh size", fontsize=25)
#leg = plt.legend(loc='best', fontsize=25)


ax.loglog(matrix[:,0], matrix[:,1], 'k', linewidth=2)
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#plt.xticks(gridSize)
ax.set_xlabel(r"$\alpha = \lambda \sigma^2 / 2$", fontsize=25)
ax.set_ylabel(r"$J = 0.5 \sum_i (T[g](\mathbf{x}_i) - \hat{T}(\mathbf{x}_i))^2$", fontsize=25)
plt.grid()
plt.show()

