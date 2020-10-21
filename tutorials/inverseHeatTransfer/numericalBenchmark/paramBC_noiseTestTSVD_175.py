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

relError_L2norm = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCnoiseTest_TSVD_175/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCnoiseTest_TSVD_175/relError_LinfNorm_mat.txt")
TSVDtruc = np.loadtxt("coarseMesh/ITHACAoutput/parameterizedBCnoiseTest_TSVD_175/TSVDtrucVector_mat.txt")

relErr_L2norm = np.empty([len(TSVDtruc)])
relErr_L2normMin = np.empty([len(TSVDtruc)])
relErr_L2normMax = np.empty([len(TSVDtruc)])
relErr_LinfNorm = np.empty([len(TSVDtruc)])
relErr_LinfNormMin = np.empty([len(TSVDtruc)])
relErr_LinfNormMax = np.empty([len(TSVDtruc)])

for i in range(len(TSVDtruc)):
    vec = relError_L2norm[i::len(TSVDtruc)]
    relErr_L2norm[i] = vec.mean() 
    relErr_L2normMin[i] = np.quantile(vec, 0.1) 
    relErr_L2normMax[i] = np.quantile(vec, 0.9) 
    vec = relError_LinfNorm[i::len(TSVDtruc)]
    relErr_LinfNorm[i] = vec.mean() 
    relErr_LinfNormMin[i] = np.quantile(vec, 0.1) 
    relErr_LinfNormMax[i] = np.quantile(vec, 0.9) 

print relErr_L2norm
f = plt.figure(2,figsize=(12,8))
plt.errorbar(TSVDtruc, relErr_L2norm, yerr=[relErr_L2norm - relErr_L2normMin, relErr_L2normMax - relErr_L2norm], markersize=15,fmt='bo', capsize=12, label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.errorbar(TSVDtruc, relErr_LinfNorm, yerr=[relErr_LinfNorm - relErr_LinfNormMin, relErr_LinfNormMax - relErr_LinfNorm], markersize=15,fmt='kv', capsize=12, label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel(r'$\alpha_{TSVD}$', fontsize=25)
plt.title(r'Noise standard dev. = $0.175$', fontsize=25)
plt.yscale('log')
plt.xlim(0,TSVDtruc[-1])
plt.grid()
plt.legend(fontsize=25)
plt.show()
