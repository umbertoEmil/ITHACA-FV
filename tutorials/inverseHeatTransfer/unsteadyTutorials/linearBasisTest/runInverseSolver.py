# runInverseSolver.py>

import files
import os

# Function that, given the mesh, the timestep size and a value of the cost function parameter
# returns the mean of the measurements discrepancy
def runInverseSolver(mesh, Deltat, costFunctionParameter, cleanCase):
    if cleanCase == 1:
        os.system("./cleanCase")

    files.sed_variable("deltaT","./system/controlDict", str(Deltat))
    files.sed_count("deltaX", mesh, "./system/blockMeshDict_input", "./system/blockMeshDict")
    files.sed_variable("costFunctionParameter","./system/ITHACAdict", str(costFunctionParameter[0]))

    # Change if you want to do a noise test
    files.sed_variable("regularizationTechnique","./constant/regularizationDict", "fullPivLU")
    files.sed_variable("addNoise","./system/ITHACAdict", "0")

    os.system("blockMesh")

    os.system("linearBasisTest")

    test_file = open('./ITHACAoutput/testInverse/measurementsDiscrepancy_mean_mat.txt', 'r')
    test_lines = test_file.readlines()
    test_file.close()
    ferNum = test_lines[0]
    return (float(ferNum))
