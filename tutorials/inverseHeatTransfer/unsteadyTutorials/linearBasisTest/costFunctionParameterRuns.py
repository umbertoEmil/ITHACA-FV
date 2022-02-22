import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import shutil

folder = "./resultsForPlots/costFunctionParameter/"
folderIn = "./ITHACAoutput/testInverse/"

#delta_x = ["100 20 85", "50 20 45", "40 15 35", "30 10 25", "20 5 15"]
#delta_t = [0.1, 0.2, 0.25, 0.5

delta_x = "40 15 35"
delta_t = "0.2"

costParam = np.geomspace(1e-16, 1e-6, 100)
print costParam
np.savetxt(folder + "/costFunctionParameter.txt", costParam, fmt='%1.4e')

heatFluxRelErr_L2_max = []
heatFluxRelErr_Linf_max = []
measurementsDiscrepancy_max = []

heatFluxRelErr_L2_mean = []
heatFluxRelErr_Linf_mean = []
measurementsDiscrepancy_mean = []

heatFluxNorm_max = []
heatFluxNorm_mean = []

with open(folder + "deltaX.txt", "w") as file:
    file.write(delta_x)

with open(folder + "deltaT.txt", "w") as file:
    file.write(delta_t)

files.sed_variable("deltaT","./system/controlDict", delta_t)
files.sed_count("deltaX", delta_x, "./system/blockMeshDict_input", "./system/blockMeshDict")
os.system("./cleanCase")
os.system("blockMesh")

files.sed_variable("regularizationTechnique","./constant/regularizationDict", "fullPivLU")
files.sed_variable("addNoise","./system/ITHACAdict", "0")

for k in costParam:
    files.sed_variable("costFunctionParameter","./system/ITHACAdict",str(k))
    os.system("linearBasisTest")
    
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

    test_file = open('./ITHACAoutput/testInverse/measurementsDiscrepancy_max_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    measurementsDiscrepancy_max.append(float(ferNum))

    test_file = open('./ITHACAoutput/testInverse/measurementsDiscrepancy_mean_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    measurementsDiscrepancy_mean.append(float(ferNum))

    test_file = open('./ITHACAoutput/testInverse/heatFluxL2norm_max_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    heatFluxNorm_max.append(float(ferNum))

    test_file = open('./ITHACAoutput/testInverse/heatFluxL2norm_mean_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    heatFluxNorm_mean.append(float(ferNum))


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

with open(folder + "measurementsDiscrepancy_max_list.txt", "w") as file:
    file_lines = "\n".join(map(str, measurementsDiscrepancy_max))
    file.write(file_lines)

with open(folder + "measurementsDiscrepancy_mean_list.txt", "w") as file:
    file_lines = "\n".join(map(str, measurementsDiscrepancy_mean))
    file.write(file_lines)

with open(folder + "heatFluxL2norm_max_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxNorm_max))
    file.write(file_lines)

with open(folder + "heatFluxL2norm_mean_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxNorm_mean))
    file.write(file_lines)
