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
    constantBasisTest_ROM.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
//#include "inverseLaplacianProblem.H"
#include "inverseHeatTransferProblem.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "ReducedSequentialIHTP_incrementalPOD_constant.H"
#include "sequentialIHTP.H"
#include "constantBasisTest_ROM.H"
#include "sequentialAcquisitionTest_steady.H"


using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    constantBasisTest_ROM example(argc, argv);
    sequentialAcquisitionTest_steady exampleSteady(argc, argv);


    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    scalar thermalCond  = 
        para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    scalar density = para->ITHACAdict->lookupOrDefault<double>("density", 0);
    scalar specificHeat = para->ITHACAdict->lookupOrDefault<double>("specificHeat",
                          0);
    scalar diffusivity = thermalCond / (density * specificHeat);
    M_Assert( diffusivity > 0, "diffusivity not specified");
    example.setDiffusivity(diffusivity);
    example.thermalCond = thermalCond;
    exampleSteady.k = thermalCond;
    M_Assert(example.thermalCond > 0, "thermalConductivity, k, not specified");
    example.density = density;
    M_Assert(example.density > 0, "Density not specified");
    example.specificHeat = specificHeat;
    M_Assert(example.specificHeat > 0, "specificHeat not specified");
    example.HTC = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    exampleSteady.H = example.HTC;
    M_Assert(example.HTC > 0, "Heat transfer coeff, H, not specified");
    example.a = para->ITHACAdict->lookupOrDefault<scalar>("a", 0);
    example.b = para->ITHACAdict->lookupOrDefault<scalar>("b", 0);
    example.c = para->ITHACAdict->lookupOrDefault<scalar>("c", 0);
    example.d = para->ITHACAdict->lookupOrDefault<scalar>("d", 0);

    example.linearTrueHeatFlux = 
        para->ITHACAdict->lookupOrDefault<bool>("linearTrueHeatFlux", 0);
    scalar timeGrad =
        para->ITHACAdict->lookupOrDefault<scalar>("timeGrad", 0);

    example.NmodesT_ic = para->ITHACAdict->lookupOrDefault<int>("NmodesT_ic", 0);
    example.T_ic_projectionTol = 
        para->ITHACAdict->lookupOrDefault<double>("T_ic_projectionTol", 0);
    int NmagicPoints =
        para->ITHACAdict->lookupOrDefault<int>("NmagicPoints", 0);
    scalar shapeParameter =
        para->ITHACAdict->lookupOrDefault<scalar>("shapeParameter", 1);
    example.maxFrequency = para->ITHACAdict->lookupOrDefault<scalar>("maxFrequency",
                           0);
    double SVDtol =
        para->ITHACAdict->lookupOrDefault<double>("SVDtol", 0);
    word PODnorm =
        para->ITHACAdict->lookupOrDefault<word>("PODnorm", "L2");

    unsigned FOMtest =
        para->ITHACAdict->lookupOrDefault<unsigned>("FOMtest", 0);
    unsigned ROMtest =
        para->ITHACAdict->lookupOrDefault<unsigned>("ROMtest", 0);

    unsigned addNoise = 
        para->ITHACAdict->lookupOrDefault<unsigned>("addNoise", 0);
    double noiseLevel =
        para->ITHACAdict->lookupOrDefault<double>("noiseLevel", 0);
    Info << "\n ************************************************************ \n";
    Info << "Conducting chirp test to compare performance of steady and unsteady inverse solvers\n";
    Info << "We assume the heat flux to estimate has the shape:\n";
    Info << "\n g = A + B sin [2 pi (a t) t] \n";
    Info << "Maximum frequency is " << example.maxFrequency << "Hz\n \n";

    example.readThermocouples();
    example.set_trueHeatFlux(timeGrad);

    /// Making a steady run for t=0 to have the unsteady initial field
    example.set_Tf(0);
    exampleSteady.restart();
    exampleSteady.set_Tf(example.Tf);
    exampleSteady.g = example.trueHeatFlux[0];
    exampleSteady.solveTrue();
    volScalarField& T(exampleSteady._T());

    for (label i = 0; i < T.internalField().size(); i++)
    {
        example.initialField.ref()[i] = T.internalField()[i];
    }

    example.solveTrue();

    word outputFolderFULL = "./ITHACAoutput/testInverse/FOM/";
    example.assignTrueIF();
    example.setParametrizedHeatFlux("rbf", shapeParameter);
    example.offlinePhase();
    
    //Set T_ic
    volScalarField initialField = example.Ttrue[0];
    
    if(addNoise)
    {
        example.addNoise(noiseLevel);
    }

    // Full 
    if(FOMtest)
    {
        example.sequentialIHTP_constant::computeHeatFluxWeights(initialField, 
                example.projectHeatFlux(example.trueHeatFlux[0]), outputFolderFULL); 
        example.inverseProblemPostProcess(outputFolderFULL);
    }

    // Reduced 
    if(ROMtest)
    {
        example.T_ic_field.resize(0);
        word outputFolderRED = "./ITHACAoutput/testInverse/ROM/";
        example.computeHeatFluxWeights(initialField,
                example.projectHeatFlux(example.trueHeatFlux[0]), outputFolderRED,
                NmagicPoints, SVDtol, PODnorm); 
        example.inverseProblemPostProcess(outputFolderRED);
        example.inverseProblemPostProcess(outputFolderRED, outputFolderFULL);
    }

    return 0;
}

