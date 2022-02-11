import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools

delta_t = [0.1, 0.2, 0.25, 0.5]


#modes_T = [1,2,3,4,5,10,15,20,25,30,35,40]
#modes_DEIM = [1,2,3,4,5,10,15,20,25,30,35,40]
#a = list(itertools.product(modes_T, modes_DEIM))
heatFluxRelErr_L2 = []
heatFluxRelErr_Linf = []

for k in delta_t:
     files.sed_variable("deltaT","./system/controlDict",str(k))
     os.system("rm -r ITHACAoutput")
     os.system("constantBasisTest")

     test_file = open('./heatFluxRelErr_L2_mat.txt', 'r')
     test_lines = test_file.readlines()
     test_file.close()
     ferNum= test_lines[0]
     heatFluxRelErr_L2.append(float(ferNum))

     test_file = open('./heatFluxRelErr_Linf_mat.txt', 'r')
     test_lines = test_file.readlines()
     test_file.close()
     ferNum= test_lines[0]
     heatFluxRelErr_Linf.append(float(ferNum))

with open("heatFluxRelErr_L2_list.txt", "w") as file:
    values = '\n'.join(map(str, heatFluxRelErr_L2))
    file.write(values)

with open("heatFluxRelErr_Linf_list.txt", "w") as file:
    file_lines = "\n".join(map(str, heatFluxRelErr_Linf))
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
