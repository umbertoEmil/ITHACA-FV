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

N = 10000
t0 = 0
tF = 10
t = np.linspace(t0, tF, N, endpoint=True)

fmin = 0.5
G = 1
f1 = fmin / tF * t
g1 = G * np.sin(2 * np.pi * f1 * t) 

f2 = 0.1 * t
g2 = G * np.sin(2 * np.pi * f2 * t) 

#with open("./ITHACAoutput/directProblem/dataAtHotSideCenter.csv") as f:
#    lis = [line.split() for line in f]        # create a list of lists
#    for i, x in enumerate(lis):              #print the list items 
#        print "line{0} = {1}".format(i, x)

relErrL2norm_unsteady = np.loadtxt("./ITHACAoutput/testInverse/relErrL2norm_mat.txt")
relErrL2norm_steady = np.loadtxt("./ITHACAoutput/steadyTest/relError_L2norm_mat.txt")

Nmeas = relErrL2norm_steady.size
t0meas = 0.
tMeas = np.linspace(t0meas, tF, Nmeas, endpoint=True)

Nunsteady = relErrL2norm_unsteady.size 
t0 = 0
tF = 10
tUnsteady = np.linspace(t0, tF, Nunsteady, endpoint=True)

fig = plt.figure(1,figsize=(8,6))
plt.plot(tUnsteady, relErrL2norm_unsteady, 'k', linewidth = 2, label='Unsteady')
plt.plot(tMeas, relErrL2norm_steady, 'ob', label='Steady')

#f = plt.figure(1,figsize=(8,6))
#plt.plot(t, g1)
#plt.plot(t, f1)
plt.xlabel('Time [s]', fontsize=25)
plt.ylabel('L2 norm of relative error', fontsize=25)
plt.grid()
plt.yscale('log')
##plt.xticks(iteration)
##plt.ylim(1e-7, 1)
##plt.xticks(np.rint(np.linspace(iteration[0], iteration[-1], num=8)))
#plt.grid()
##plt.savefig('J.eps', bbox_inches='tight')
plt.legend()

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (s)', fontsize=25)
ax1.set_ylabel('Frequency [Hz]', color=color, fontsize=25)  # we already handled the x-label with ax1
ax1.plot(t, f1, color=color, linewidth=2)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Amplitude', color=color, fontsize=25)
ax2.plot(t, g1, color=color, linewidth=2)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

