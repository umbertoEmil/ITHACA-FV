import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import shutil

delta_x = ["100 20 85", "50 20 45", "40 15 35", "30 10 25", "20 5 15"]
delta_t = [0.1, 0.2, 0.25, 0.5]
a = list(itertools.product(delta_x, delta_t))
print(a)

heatFluxRelErr_L2_max = []
heatFluxRelErr_Linf_max = []

heatFluxRelErr_L2_mean = []
heatFluxRelErr_Linf_mean = []

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


for k in a:
    files.sed_count("deltaX", str(k[0]), "./system/blockMeshDict_input", "./system/blockMeshDict")
    files.sed_variable("deltaT","./system/controlDict",str(k[1]))
    os.system("./cleanCase")
    os.system("blockMesh")
    os.system("constantBasisTest")
    shutil.copyfile("./ITHACAoutput/testInverse/heatFluxRelErr_L2_mat.txt", "./resultsForPlots/heatFluxRelErr_L2_" + str(k[0]) + "_" + str(k[1]))
    shutil.copyfile("./ITHACAoutput/testInverse/heatFluxRelErr_Linf_mat.txt", "./resultsForPlots/heatFluxRelErr_Linf_" + str(k[0]) + "_" + str(k[1]) )
    shutil.copyfile("./ITHACAoutput/testInverse/timeVec_mat.txt", "./resultsForPlots/timeVec_" + str(k[1]) )
    shutil.copyfile("./ITHACAoutput/testInverse/probeHeat_true_mat.txt", "./resultsForPlots/probeHeat_true_" + str(k[1]) )
    shutil.copyfile("./ITHACAoutput/testInverse/probeHeat_rec_mat.txt", "./resultsForPlots/probeHeat_rec_" + str(k[0]) + "_" + str(k[1]) )
        
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
    
    test_file = open('./onlineCPUtime_avg_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    CPUtime.append(float(ferNum))

with open("./resultsForPlots/heatFluxRelErr_L2_max_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2_max))
    file.write(values)

with open("./resultsForPlots/heatFluxRelErr_Linf_max_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf_max))
    file.write(file_lines)

with open("./resultsForPlots/heatFluxRelErr_L2_mean_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2_mean))
    file.write(values)

with open("./resultsForPlots/heatFluxRelErr_Linf_mean_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf_mean))
    file.write(file_lines)

with open("./resultsForPlots/CPUtime.txt", "w") as file:
    file_lines = "\n".join(map(str, CPUtime))
    file.write(file_lines)

