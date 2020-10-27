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

relError_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/ParamBCnoiseLevelTest/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/ParamBCnoiseLevelTest/relError_LinfNorm_mat.txt")

#Have a look at the vector noiseLevel in analyticalBenchmark.C 
noiseLevel = np.array([.005, .01, .02, .03, .04, .05, .06, .07, .08, .1])

Ntests = int(len(relError_L2norm) / len(noiseLevel))

relErr_L2norm = np.empty([len(noiseLevel)])
relErr_LinfNorm = np.empty([len(noiseLevel)])
for i in range(len(noiseLevel)):
    vec = relError_L2norm[i * Ntests:(i+1) * Ntests - 1]
    relErr_L2norm[i] = vec.mean() 
    vec = relError_LinfNorm[i * Ntests:(i+1) * Ntests - 1]
    relErr_LinfNorm[i] = vec.mean() 

print relErr_L2norm
f = plt.figure(3,figsize=(12,8))
plt.semilogy(noiseLevel, relErr_L2norm, 'bo', markersize=15, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.semilogy(noiseLevel, relErr_LinfNorm, 'kv', markersize=15, label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel('Noise standard deviation', fontsize=25)
plt.ylabel('Mean of relative error norms', fontsize=25)
plt.grid()
plt.title(r"LU", fontsize=25)
plt.legend(fontsize=25)
plt.show()
