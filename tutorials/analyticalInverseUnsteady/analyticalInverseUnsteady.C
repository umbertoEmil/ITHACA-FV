/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    Example of a heat transfer Reduction Problem
SourceFiles
    DanieliMold.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseHeatTransferProblem.H"
// #include "reducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
//#include "ITHACAbayesian.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "sampledTriSurfaceMesh.H"

#include "analyticalInverseUnsteady.H"
using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    // Reading tests to perform
    ITHACAparameters para;
    scalar k = para.ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    scalar density = para.ITHACAdict->lookupOrDefault<double>("density", 0);
    scalar specificHeat = para.ITHACAdict->lookupOrDefault<double>("specificHeat", 0);
    scalar diffusivity = k / (density * specificHeat);
    
    analyticalInverseUnsteady example(argc, argv, diffusivity);
    
    M_Assert( diffusivity > 0, "diffusivity not specified");
    example.k = k;
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.density = density;
    M_Assert(example.density > 0, "Density not specified");
    example.specificHeat = specificHeat;
    M_Assert(example.specificHeat > 0, "specificHeat not specified");
    example.H = para.ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    example.exponentialSolution = para.ITHACAdict->lookupOrDefault<unsigned>("exponentialSolution", 0);
    example.quadraticSolution = para.ITHACAdict->lookupOrDefault<unsigned>("quadraticSolution", 0);
    M_Assert(example.exponentialSolution || example.quadraticSolution,
        "Chose one between quadraticSolution and exponentialSolution");

    unsigned directProblemTest = para.ITHACAdict->lookupOrDefault<unsigned>("directProblemTest", 0);
    unsigned parameterizedBCtest = para.ITHACAdict->lookupOrDefault<unsigned>("parameterizedBCtest", 0);
    
    
    //example.heatFlux = cnpy::load(example.heatFlux, "DanieliHeatFluxAndCastingSpeed.npy");
    //std::cout << "heatFlux.rows = " << example.heatFlux.rows() << std::endl;

    example.readThermocouples();
    example.analyticalSolution();
    if(directProblemTest)
    {
        example.solveTrue();
        example.directProblemPostProcess();
    }

    if(parameterizedBCtest)
    {
        word outputFolder = "./ITHACAoutput/testInverse/";
        example.assignTrueIF();
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        example.parameterizedBC(outputFolder, "fullPivLU");
	example.inverseProblemPostProcess(outputFolder);
    }
    //example.reconstrucT();
    //example.solveUnsteady();
    
    return 0;
}

