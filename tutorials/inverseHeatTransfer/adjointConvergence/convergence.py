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

from numpy import genfromtxt
my_data = genfromtxt('relErr.csv', delimiter=',')
cellSize = np.array([0.25, 0.1, 0.0667, 0.03333, 0.025])

relErr = my_data[:,1]
print my_data

C0 = relErr[-1] / ( cellSize[-1] * cellSize[-1] )
quadraticLine = C0  *  np.power(cellSize,2)
C0 = relErr[-1] / ( cellSize[-1] )
linearLine = C0  *  cellSize

#C1 = relDifferenceL2[-1] / ( cellSize[-1] * cellSize[-1] )
#relQuadraticLine = C1 * np.power(cellSize,2)

f = plt.figure(2,figsize=(13,9))
ax1 = plt.axes(xscale='log', yscale='log')
ax1.plot(cellSize,  my_data[:,1], 'bo', markersize=15, label=r'Relative error')
ax1.plot(cellSize, quadraticLine, 'b--', linewidth=2, label=r'Quadratic')
ax1.plot(cellSize, linearLine, 'b-.', linewidth=2, label=r'Linear')
#
#ax1.plot(cellSize, relDifferenceL2, 'kv', markersize=15, label=r'Relative error')
#ax1.plot(cellSize, relQuadraticLine, 'k-', linewidth=2, label="Quadratic convergence")

ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
plt.xticks(cellSize)
plt.xlabel("Cell edge length [m]", fontsize=25)
plt.legend(loc='best', fontsize=25)
plt.grid()


plt.show()
