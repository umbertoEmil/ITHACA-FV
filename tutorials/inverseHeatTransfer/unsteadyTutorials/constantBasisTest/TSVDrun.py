import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import files
import itertools
import shutil

alpha = np.linspace(5, 100, 95, dtype = int, endpoint=True)
print(alpha)

heatFluxRelErr_L2 = []
heatFluxRelErr_Linf = []

with open("./resultsForPlots/TSVD/TSVDalpha.txt", "w") as file:
    values = '\n'.join(map(str, alpha)) 
    file.write(values)

os.system("./cleanCase")
os.system("blockMesh")
files.sed_variable("regularizationTechnique", "./constant/regularizationDict", "TSVD")

for k in alpha:
     files.sed_variable("TSVD_filter", "./constant/regularizationDict", k)
     os.system("constantBasisTest")
     shutil.copyfile("./ITHACAoutput/testInverse/heatFluxRelErr_L2_mat.txt", "./resultsForPlots/TSVD/heatFluxRelErr_L2_TSVD_" + str(k))
     shutil.copyfile("./ITHACAoutput/testInverse/heatFluxRelErr_Linf_mat.txt", "./resultsForPlots/TSVD/heatFluxRelErr_Linf_TSVD_" + str(k))
     shutil.copyfile("./ITHACAoutput/testInverse/timeVec_mat.txt", "./resultsForPlots/TSVD/timeVec")
     shutil.copyfile("./ITHACAoutput/testInverse/probeHeat_true_mat.txt", "./resultsForPlots/TSVD/probeHeat_true")
     shutil.copyfile("./ITHACAoutput/testInverse/probeHeat_rec_mat.txt", "./resultsForPlots/TSVD/probeHeat_rec_TSVD_" + str(k))

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

with open("./resultsForPlots/TSVD/heatFluxRelErr_L2_TSVD.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2))
    file.write(values)

with open("./resultsForPlots/TSVD/heatFluxRelErr_Linf_TSVD.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf))
    file.write(file_lines)
