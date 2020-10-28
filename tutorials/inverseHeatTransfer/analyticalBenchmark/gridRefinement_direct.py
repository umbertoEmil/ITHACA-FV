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


differenceL2 = np.array([
    0.0359046014773,
    0.00896508653386,
    0.0039724735755,
    0.00222477783325,
    0.00097633449241,
    0.000539361747386
])

relDifferenceL2 = np.array([
    0.00138740460926,
    0.000346423630056,
    0.000153502000358,
    8.59685637337e-05,
    3.77269463861e-05,
    2.08417011634e-05,
])


f = plt.figure(2,figsize=(12,8))
ax1 = plt.axes(xscale='log', yscale='log')
ax1.plot(cellSize, differenceL2, 'bo', markersize=15, label=r'Error')
ax1.plot(cellSize, relDifferenceL2, 'kv', markersize=15, label=r'Relative error')
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
plt.xticks(cellSize)
plt.xlabel("Cell edge lenght [m]", fontsize=25)
plt.legend(loc='best', fontsize=25)
plt.grid()



plt.show()
