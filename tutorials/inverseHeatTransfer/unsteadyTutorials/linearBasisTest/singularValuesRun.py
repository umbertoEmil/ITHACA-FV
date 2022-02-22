import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import shutil

delta_x = ["100 20 85", "50 20 45", "40 15 35", "30 10 25", "20 5 15"]
#delta_x = ["160 8 96", "100 5 60", "80 4 42"]
delta_t = [0.1, 0.2, 0.25, 0.5]
a = list(itertools.product(delta_x, delta_t))
print(a)
heatFluxRelErr_L2 = []
heatFluxRelErr_Linf = []
CPUtime = []

with open("./resultsForPlots/deltaX.txt", "w") as file:
    values = '\n'.join(delta_x)
    file.write(values)

with open("./resultsForPlots/deltaT.txt", "w") as file:
    values = '\n'.join(map(str, delta_t)) 
    file.write(values)

with open("./resultsForPlots/deltaXdeltaT.txt", "w") as file:
    values = '\n'.join(map(str, a)) 
    file.write(values)

files.sed_variable("conditioningTest","./system/ITHACAdict",str(1))
files.sed_variable("reconstructionTest","./system/ITHACAdict",str(0))
files.sed_variable("thermocoupleReconstructionTest","./system/ITHACAdict",str(0))
files.sed_variable("inverseTest","./system/ITHACAdict",str(0))

for k in a:
     files.sed_count("deltaX", str(k[0]), "./system/blockMeshDict_input", "./system/blockMeshDict")
     files.sed_variable("deltaT","./system/controlDict",str(k[1]))
     os.system("./cleanCase")
     os.system("blockMesh")
     os.system("linearBasisTest")
     shutil.copyfile("./ITHACAoutput/conditioningTest/singularValues_mat.txt", "./resultsForPlots/conditioningTest/singularValues_" + str(k[0]) + "_" + str(k[1]))
     shutil.copyfile("./ITHACAoutput/conditioningTest/conditioningNumber_mat.txt", "./resultsForPlots/conditioningTest/conditioningNumber_" + str(k[0]) + "_" + str(k[1]))

