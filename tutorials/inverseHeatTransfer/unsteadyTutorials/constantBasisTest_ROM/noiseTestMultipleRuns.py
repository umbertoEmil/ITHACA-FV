import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import shutil


####################################
####### USER DEFINED INPUTS ########
####################################
solver = "fullPivLU"
TSVD_filter = 10
delta_x = "40 15 35"
delta_t = "0.5"
Ntests = 1
####################################

heatFluxRelErr_L2 = []
heatFluxRelErr_Linf = []

heatFluxRelErr_L2_max = []
heatFluxRelErr_Linf_max = []

heatFluxRelErr_L2_mean = []
heatFluxRelErr_Linf_mean = []

folder = "./resultsForPlots/noiseTest/" + solver + "/"
folderIn = "./ITHACAoutput/testInverse/ROM/"


with open(folder + "deltaX.txt", "w") as file:
    file.write(delta_x)

with open(folder + "deltaT.txt", "w") as file:
    file.write(delta_t)

with open(folder + "solver.txt", "w") as file:
    file.write(solver)

with open(folder + "TSVD_filter.txt", "w") as file:
    file.write(str(TSVD_filter))

files.sed_variable("deltaT","./system/controlDict", delta_t)
files.sed_count("deltaX", delta_x, "./system/blockMeshDict_input", "./system/blockMeshDict")
os.system("./cleanCase")
os.system("blockMesh")

files.sed_variable("regularizationTechnique","./constant/regularizationDict", solver)
files.sed_variable("TSVD_filter","./constant/regularizationDict", TSVD_filter)
files.sed_variable("addNoise","./system/ITHACAdict", "1")
files.sed_variable("costFunctionParameter","./system/ITHACAdict", "0.0")
files.sed_variable("linearTrueHeatFlux","./system/ITHACAdict", "0")

noiseLevel = np.geomspace(1e-4, 1e-2, 2)


with open(folder + "noiseLevel.txt", "w") as file:
    values = '\n'.join(map(str, noiseLevel))
    file.write(values)



for k in noiseLevel:
    files.sed_variable("noiseStdDev","./system/ITHACAdict",str(k))
    for n in range(Ntests):

        print("**********************************************")
        print("\n\nTest " + str(n) + ".\n\nnoiseLevel = " + str(k) + "\n\n")
        print("**********************************************")
        os.system("constantBasisTest_ROM")
    
        test_file = open('./heatFluxRelErr_L2_max_mat.txt', 'r')
        test_lines = test_file.readlines()
        test_file.close()
        ferNum = test_lines[0]
        heatFluxRelErr_L2_max.append(float(ferNum))

        test_file = open('./heatFluxRelErr_Linf_max_mat.txt', 'r')
        test_lines = test_file.readlines()
        test_file.close()
        ferNum = test_lines[0]
        heatFluxRelErr_Linf_max.append(float(ferNum))
        
        test_file = open('./heatFluxRelErr_L2_mean_mat.txt', 'r')
        test_lines = test_file.readlines()
        test_file.close()
        ferNum = test_lines[0]
        heatFluxRelErr_L2_mean.append(float(ferNum))

        test_file = open('./heatFluxRelErr_Linf_mean_mat.txt', 'r')
        test_lines = test_file.readlines()
        test_file.close()
        ferNum = test_lines[0]
        heatFluxRelErr_Linf_mean.append(float(ferNum))

with open(folder + "heatFluxRelErr_L2_max_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2_max))
    file.write(values)

with open(folder + "heatFluxRelErr_Linf_max_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf_max))
    file.write(file_lines)

with open(folder + "heatFluxRelErr_L2_mean_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2_mean))
    file.write(values)

with open(folder + "heatFluxRelErr_Linf_mean_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf_mean))
    file.write(file_lines)

