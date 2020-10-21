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
#deltaT = [0.001, 0.003, 0.005, 0.008, 0.01, 0.03, 0.05, 0.08, 0.1]        
deltaT = [0.001, 0.01, 0.1]
Ncells = [20]
diffusivity = 2.0 / (3.0 * 2000.0)
print"Diffusivity = "
print(diffusivity)

relErrL2 = np.zeros((len(Ncells), len(deltaT)))
relErrLinf = np.zeros((len(Ncells), len(deltaT)))
CFL = np.zeros((len(Ncells), len(deltaT)))
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
         Dx = 1. / meshI
         courantNo = diffusivity * k * 3. / ( Dx * Dx )
         CFL.itemset((i,j), courantNo)
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
np.save("./CFL", CFL)
