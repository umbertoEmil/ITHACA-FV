import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
import numpy as np 
import sys
from scipy import linalg
from scipy.stats.kde import gaussian_kde
from numpy import linspace
import seaborn as sns


#plt.style.use('classic')
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

seeds = np.loadtxt  ("./seeds_mat.txt") #L2 norm of relative error
pdf = np.loadtxt  ("./pdf_mat.txt") #L2 norm of relative error
out = np.loadtxt  ("./output_mat.txt") #L2 norm of relative error

f = plt.figure(20,figsize=(12,8))
        sns.distplot(out, norm_hist = True, kde=True,
                     color = 'darkblue',
                     hist_kws={'edgecolor':'black'},
                     kde_kws={'linewidth': 4})

#plt.plot( out, 'o', linewidth=2)

plt.grid()

#g = plt.figure(10,figsize=(12,8))
#plt.plot( seeds, pdf, 'o', linewidth=2)
#plt.ylabel(r"pdf", fontsize=25)
#plt.xlabel(r"seeds", fontsize=25)
#plt.grid()
plt.show()
