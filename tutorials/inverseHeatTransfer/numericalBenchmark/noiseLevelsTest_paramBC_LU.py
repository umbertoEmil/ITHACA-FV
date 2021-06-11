import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relError_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/ParamBCnoiseLevelTest_LU/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/ParamBCnoiseLevelTest_LU/relError_LinfNorm_mat.txt")

noiseLevel = np.loadtxt("coarseMesh/ITHACAoutput/ParamBCnoiseLevelTest_LU/noiseLevelVec_mat.txt")

print noiseLevel

Ntests = int(len(relError_L2norm) / len(noiseLevel))

relErr_L2norm = np.empty([len(noiseLevel)])
relErr_L2normMin = np.empty([len(noiseLevel)])
relErr_L2normMax = np.empty([len(noiseLevel)])
relErr_LinfNorm = np.empty([len(noiseLevel)])
relErr_LinfNormMin = np.empty([len(noiseLevel)])
relErr_LinfNormMax = np.empty([len(noiseLevel)])
for i in range(len(noiseLevel)):
    vec = relError_L2norm[i * Ntests:(i+1) * Ntests - 1]
    relErr_L2norm[i] = vec.mean() 
    relErr_L2normMin[i] = np.quantile(vec, 0.1) 
    relErr_L2normMax[i] = np.quantile(vec, 0.9) 
    vec = relError_LinfNorm[i * Ntests:(i+1) * Ntests - 1]
    relErr_LinfNorm[i] = vec.mean() 
    relErr_LinfNormMin[i] = np.quantile(vec, 0.1) 
    relErr_LinfNormMax[i] = np.quantile(vec, 0.9) 

print relErr_L2norm
#f = plt.figure(3,figsize=(12,8))
#plt.semilogy(noiseLevel, relErr_L2norm, 'bo', markersize=15, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
#plt.semilogy(noiseLevel, relErr_LinfNorm, 'kv', markersize=15, label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
#plt.xlabel(r'Noise standard deviation on $G$', fontsize=25)
#plt.ylabel('Mean of relative error norms', fontsize=25)
#plt.grid(True, which="both", ls="-")
#plt.title(r"LU", fontsize=25)
#plt.legend(fontsize=25)

f = plt.figure(5,figsize=(12,8))
plt.errorbar(noiseLevel, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normMax - relErr_L2norm], markersize=15,fmt='bo', capsize=12, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.errorbar(noiseLevel, relErr_LinfNorm, yerr=[relErr_LinfNorm - relErr_LinfNormMin, relErr_LinfNormMax - relErr_LinfNorm], markersize=15,fmt='kv', capsize=12, label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel(r'Noise standard deviation', fontsize=25)
plt.ylabel('Mean of relative error norms', fontsize=25)
plt.grid(True, which="both", ls="-")
#plt.title(r"LU", fontsize=25)
#plt.title(r"TSVD, $\alpha_{TSVD} = 7$", fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0)
plt.legend(fontsize=25)
ax = plt.gca()
ax.set_yscale('log')
#plt.tick_params(axis='y', which='minor')
#ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0e"))

plt.show()
