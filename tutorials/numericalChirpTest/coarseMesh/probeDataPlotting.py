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

N = 101
t0 = 0
tF = 10
t = np.linspace(t0, tF, N, endpoint=True)

fmin = 0.5
G = 1
f1 = fmin / tF * t
g1 = G * np.sin(2 * np.pi * f1 * t) 

from numpy import genfromtxt
my_data = genfromtxt('/media/umberto/OS/Users/Umberto/Documents/PhD/ITHACA-FV/tutorials/numericalChirpTest/coarseMesh/ITHACAoutput/directProblem/probeData0.csv', delimiter=',')

T = my_data[1:, 2]

my_data = genfromtxt('/media/umberto/OS/Users/Umberto/Documents/PhD/ITHACA-FV/tutorials/numericalChirpTest/coarseMesh/ITHACAoutput/directProblem/probeDataAtFace0.csv', delimiter=',')

gComp = my_data[1:, 9]

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (s)', fontsize=25)
ax1.set_ylabel('Temperature at thermocouple [K]', color=color, fontsize=25)  # we already handled the x-label with ax1
ax1.plot(t, T, color=color, linewidth=2)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('g [W/m2]', color=color, fontsize=25)
ax2.plot(t, gComp, color=color, linewidth=2)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()

