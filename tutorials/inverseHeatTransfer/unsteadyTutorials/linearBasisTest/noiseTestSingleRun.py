import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import shutil

delta_x = "40 15 35"
delta_t = "0.5"

heatFluxRelErr_L2 = []
heatFluxRelErr_Linf = []

folder = "./resultsForPlots/noiseTest/LU/"
folderIn = "./ITHACAoutput/testInverse/"


with open(folder + "deltaX.txt", "w") as file:
    file.write(delta_x)

with open(folder + "deltaT.txt", "w") as file:
    file.write(delta_t)

files.sed_variable("deltaT","./system/controlDict", delta_t)
files.sed_count("deltaX", delta_x, "./system/blockMeshDict_input", "./system/blockMeshDict")
os.system("./cleanCase")
os.system("blockMesh")

files.sed_variable("regularizationTechnique","./constant/regularizationDict", "fullPivLU")
files.sed_variable("addNoise","./system/ITHACAdict", "1")

noiseLevel = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1]


with open(folder + "noiseLevel.txt", "w") as file:
    values = '\n'.join(map(str, noiseLevel))
    file.write(values)



for k in noiseLevel:
     files.sed_variable("noiseStdDev","./system/ITHACAdict",str(k))
     os.system("linearBasisTest")
     shutil.copyfile(folderIn + "/heatFluxRelErr_L2_mat.txt",   folder + "heatFluxRelErr_L2_noise" + str(k))
     shutil.copyfile(folderIn + "/heatFluxRelErr_Linf_mat.txt", folder + "heatFluxRelErr_Linf_noise" + str(k))
     shutil.copyfile(folderIn + "/timeVec_mat.txt",             folder + "timeVec")
     shutil.copyfile(folderIn + "/probeHeat_true_mat.txt",      folder + "probeHeat_true_noise" + str(k))
     shutil.copyfile(folderIn + "/probeHeat_rec_mat.txt",       folder + "probeHeat_rec_noise" + str(k))

     test_file = open('./heatFluxRelErr_L2_mat.txt', 'r')
     test_lines = test_file.readlines()
     test_file.close()
     ferNum = test_lines[0]
     heatFluxRelErr_L2.append(float(ferNum))

     test_file = open('./heatFluxRelErr_Linf_mat.txt', 'r')
     test_lines = test_file.readlines()
     test_file.close()
     ferNum = test_lines[0]
     heatFluxRelErr_Linf.append(float(ferNum))

with open(folder + "heatFluxRelErr_L2_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2))
    file.write(values)

with open(folder + "heatFluxRelErr_Linf_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf))
    file.write(file_lines)

