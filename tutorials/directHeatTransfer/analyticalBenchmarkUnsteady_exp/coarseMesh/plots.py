import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pylab as pylab
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

relErr_L2 = np.load("./relErrL2.npy")
relErr_Linf = np.load("./relErrLinf.npy")
Ncells = np.load("./Ncells.npy")
deltaT = np.load("./deltaT.npy")
edgeLength = 1.0 / Ncells


print deltaT
print relErr_L2
print Ncells


relErr_L2 = np.transpose(relErr_L2)
relErr_Linf = np.transpose(relErr_Linf)
fig = plt.figure(2,figsize=(12,8))
ax1 = plt.axes(xscale='log', yscale='log')
#i = 0
#for row in relErr_L2:
#    color = next(ax1._get_lines.prop_cycler)['color']
#    ax1.plot(edgeLength, row, "o", markersize=15, label=r"$\Delta t $" + '% 6.2f' % deltaT[i], color = color)
#    ax1.plot(edgeLength, relErr_Linf[i,:], "v", markersize=15, color = color)
#    i = i + 1

ax1.plot(edgeLength, relErr_L2[1,:], "bo", markersize=15, label=r"$|| \epsilon ||_{L^2(\Omega_s \times (0, t_f])}$")
ax1.plot(edgeLength, relErr_Linf[1,:], "kv", markersize=15,label=r"$|| \epsilon ||_{L^\infty(\Omega_s \times (0, t_f])}$") 

C0 = relErr_L2[0,-1] / ( edgeLength[-1] * edgeLength[-1] )
quadraticLine = C0  *  np.power(edgeLength,2)
ax1.plot(edgeLength, quadraticLine, 'b-', linewidth=2)

C0 = relErr_Linf[0,-1] / ( edgeLength[-1] * edgeLength[-1] )
quadraticLine = C0  *  np.power(edgeLength,2)
ax1.plot(edgeLength, quadraticLine, 'k-', linewidth=2, label = "Quadratic")

#ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax1.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.semilogx(gridSize,AbsError_LinfNorm, 'go')
#plt.xticks(edgeLength)
#ax1.set_xticks(edgeLength)
plt.title(r"$\Delta t = 0.1 s$", fontsize=25)
plt.xlabel(r'Cell edge length [$m$]', fontsize=25)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(loc='best', fontsize=25)

plt.show()
