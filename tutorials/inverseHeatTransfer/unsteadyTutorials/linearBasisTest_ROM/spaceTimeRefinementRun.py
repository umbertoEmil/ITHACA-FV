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


for k in a:
     files.sed_count("deltaX", str(k[0]), "./system/blockMeshDict_input", "./system/blockMeshDict")
     files.sed_variable("deltaT","./system/controlDict",str(k[1]))
     os.system("./cleanCase")
     os.system("blockMesh")
     os.system("linearBasisTest_ROM")
     shutil.copyfile("./ITHACAoutput/testInverse/ROM/heatFluxRelErr_L2_mat.txt", "./resultsForPlots/heatFluxRelErr_L2_" + str(k[0]) + "_" + str(k[1]))
     shutil.copyfile("./ITHACAoutput/testInverse/ROM/heatFluxRelErr_Linf_mat.txt", "./resultsForPlots/heatFluxRelErr_Linf_" + str(k[0]) + "_" + str(k[1]) )
     shutil.copyfile("./ITHACAoutput/testInverse/ROM/timeVec_mat.txt", "./resultsForPlots/timeVec_" + str(k[1]) )
     shutil.copyfile("./ITHACAoutput/testInverse/ROM/probeHeat_true_mat.txt", "./resultsForPlots/probeHeat_true_" + str(k[1]) )
     shutil.copyfile("./ITHACAoutput/testInverse/ROM/probeHeat_rec_mat.txt", "./resultsForPlots/probeHeat_rec_" + str(k[0]) + "_" + str(k[1]) )

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

     test_file = open('./onlineCPUtime_avg_mat.txt', 'r')
     test_lines = test_file.readlines()
     test_file.close()
     ferNum = test_lines[0]
     CPUtime.append(float(ferNum))

with open("./resultsForPlots/heatFluxRelErr_L2_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2))
    file.write(values)

with open("./resultsForPlots/heatFluxRelErr_Linf_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf))
    file.write(file_lines)

with open("./resultsForPlots/CPUtime.txt", "w") as file:
    file_lines = "\n".join(map(str, CPUtime))
    file.write(file_lines)

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
