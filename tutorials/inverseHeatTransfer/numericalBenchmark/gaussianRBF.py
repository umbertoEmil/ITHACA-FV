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
xAxis = np.linspace(-4, 4, 1000) 
eta = 1

eta = 5e-1
f = plt.figure(2,figsize=(14,8))
for i in range(5):
    rbf = np.exp( - eta * eta * xAxis * xAxis)
    eta = eta * 2 
    plt.plot(xAxis, rbf)

plt.show()
