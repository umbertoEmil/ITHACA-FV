import numpy as np
import os
import shutil
import subprocess
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools

deltaT = [0.01, 0.1, 1]
Ncells = [10, 15, 20, 25, 30, 35]

relErrL2 = np.zeros((len(Ncells), len(deltaT)))
relErrLinf = np.zeros((len(Ncells), len(deltaT)))
i = 0
for meshI in Ncells:
    shutil.copyfile('./system/blockMeshDictOrig', './system/blockMeshDict')
    files.sed_count("Ncells",str(meshI),"./system/blockMeshDict", None, 3)
    os.system("blockMesh")
    j = 0
    for k in deltaT:
         files.sed_variable("deltaT","./system/controlDict",str(k))
         os.system("rm -r ITHACAoutput")
         os.system("analyticalBenchmark_unsteady")
         L2 = np.loadtxt("./relErrL2_time.txt")
         relErrL2.itemset((i,j), L2)
         Linf = np.loadtxt("./relErrLinf_time.txt")
         relErrLinf.itemset((i,j), Linf)
         j = j + 1
    i = i + 1

print"relErrL2 = "
print(relErrL2)

print"relErrLinf = "
print(relErrLinf)

np.save("./relErrL2", relErrL2)
np.save("./relErrLinf", relErrLinf)
np.save("./deltaT", deltaT)
np.save("./Ncells", Ncells)


# error_total=[]
# error=[]
# for k,j in zip(modes_T, modes_DEIM):
#      s = "error_"+str(k)+"_"+str(j)+"_"+str(j)+"_mat.py"
#      m = "error_"+str(k)+"_"+str(j)+"_"+str(j)
#      exec(open(s).read())
#      exec("error_total.append("+m+")")
# for j in range(0,len(modes_DEIM)):
#     error.append(np.mean(error_total[j]))
# print(error)
# plt.semilogy(modes_DEIM,error,':o', label='Relative error for ROM')
# # plt.semilogy(PRO[:,0],PRO[:,1],'k--v', label='Relative error for L2 proj.')
# # plt.xlim(5,50)
# plt.xlabel("$N$ of modes")
# plt.ylabel("L2 Rel. Error.")
# # plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
# plt.grid(True)
# # f.savefig("poisson.pdf", bbox_inches='tight')
# plt.show()
