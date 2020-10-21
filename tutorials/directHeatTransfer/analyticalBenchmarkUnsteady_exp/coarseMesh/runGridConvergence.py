import numpy as np
import os
import shutil
import subprocess
import files 
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools

#deltaT = np.logspace(math.log10(0.01), math.log10(1), 8)
deltaT = [0.01, 0.1, 1]
Ncells = [10, 20, 30, 40]

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
         os.system("analyticalBenchmarkUnsteady_exp")
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
