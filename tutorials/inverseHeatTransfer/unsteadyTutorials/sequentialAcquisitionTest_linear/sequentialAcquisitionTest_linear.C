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
    sequentialAcquisitionTest_unsteady.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "inverseHeatTransferProblem.H"
#include "sequentialIHTP_linear.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "sequentialAcquisitionTestLinear_unsteady.H"
#include "sequentialAcquisitionTest_steady.H"

using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    sequentialAcquisitionTestLinear_unsteady example(argc, argv);
    sequentialAcquisitionTest_steady exampleSteady(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    scalar k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    scalar density = para->ITHACAdict->lookupOrDefault<double>("density", 0);
    scalar specificHeat = para->ITHACAdict->lookupOrDefault<double>("specificHeat",
                          0);
    scalar diffusivity = k / (density * specificHeat);
    M_Assert( diffusivity > 0, "diffusivity not specified");
    example.setDiffusivity(diffusivity);
    example.k = k;
    exampleSteady.k = k;
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.density = density;
    M_Assert(example.density > 0, "Density not specified");
    example.specificHeat = specificHeat;
    M_Assert(example.specificHeat > 0, "specificHeat not specified");
    example.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    exampleSteady.H = example.H;
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    example.a = para->ITHACAdict->lookupOrDefault<scalar>("a", 0);
    example.b = para->ITHACAdict->lookupOrDefault<scalar>("b", 0);
    example.c = para->ITHACAdict->lookupOrDefault<scalar>("c", 0);
    example.d = para->ITHACAdict->lookupOrDefault<scalar>("d", 0);
    example.T0 = para->ITHACAdict->lookupOrDefault<scalar>("T0", 0);
    example.timeBasisType = para->ITHACAdict->lookupOrDefault<word>("timeBasisType",
                            "None");
    word linSysSolver = para->ITHACAdict->lookupOrDefault<word>("linSysSolver",
                        "None");
    label TSVDtruncation =
        para->ITHACAdict->lookupOrDefault<label>("TSVDtruncation", 0);
    scalar shapeParameter =
        para->ITHACAdict->lookupOrDefault<scalar>("shapeParameter", 1);
    example.maxFrequency = para->ITHACAdict->lookupOrDefault<scalar>("maxFrequency",
                           0);
    unsigned parameterizedBC_steadyTest =
        para->ITHACAdict->lookupOrDefault<unsigned>("parameterizedBC_steadyTest", 0);
    unsigned parameterizedBC_unsteadyTest =
        para->ITHACAdict->lookupOrDefault<unsigned>("parameterizedBC_unsteadyTest", 0);
    Info << "\n ************************************************************ \n";
    Info << "Conducting chirp test to compare performance of steady and unsteady inverse solvers\n";
    Info << "We assume the heat flux to estimate has the shape:\n";
    Info << "\n g = A + B sin [2 pi (a t) t] \n";
    Info << "Maximum frequency is " << example.maxFrequency << "Hz\n \n";
    Info << "Using time basis of type " << example.timeBasisType << endl << endl;
    example.readThermocouples();
    example.set_gTrue();
    ///// Making a steady run for t=0 to have the unsteady initial field
    example.set_Tf(0);
    exampleSteady.restart();
    exampleSteady.set_Tf(example.Tf);
    exampleSteady.g = example.gTrue[0];
    exampleSteady.solveTrue();
    volScalarField& T(exampleSteady._T());

    for (label i = 0; i < T.internalField().size(); i++)
    {
        example.initialField.ref()[i] = T.internalField()[i];
    }

    example.solveTrue();

    if (parameterizedBC_steadyTest)
    {
        exampleSteady.restart();
        exampleSteady.readThermocouples();
        exampleSteady.set_Tf(example.Tf);
        exampleSteady.Tmeas.resize(example.thermocouplesNum);
        word outputFolder = "./ITHACAoutput/steadyTest/";
        Info << "*********************************************************" << endl;
        Info << "Performing offline part of STEADY test\n\n";
        exampleSteady.set_gParametrized("rbf", shapeParameter);
        exampleSteady.parameterizedBCoffline();
        Eigen::VectorXd timeVector(example.samplingTime.size());
        forAll(example.samplingTime, sampleI)
        {
            timeVector(sampleI) = example.samplingTime[sampleI];
        }
        forAll(example.samplingTime, sampleI)
        {
            Info << "*********************************************************" << endl;
            Info << "Performing STEADY test for sampling time " <<
                 example.samplingTime[sampleI] << " s\n\n";
            exampleSteady.Tmeas = example.Tmeas.segment(sampleI * example.thermocouplesNum,
                                  example.thermocouplesNum);
            scalar realTime = example.samplingTime[sampleI];
            exampleSteady.gTrue = example.gTrue[example.samplingSteps[sampleI]];
            volScalarField gTrueField = exampleSteady.list2Field(exampleSteady.gTrue);
            ITHACAstream::exportSolution(gTrueField,
                                         std::to_string(realTime), outputFolder,
                                         "gTrue");
            Eigen::VectorXd residualNorms;
            Info << "Solver: fullPivLU " << endl;
            Info << endl;
            exampleSteady.parameterizedBC("fullPivLU", 3);
            Info << "Computation ENDED, saving solution" << endl;
            Info << endl;
            volScalarField gParametrizedField = exampleSteady.list2Field(exampleSteady.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(realTime),
                                         outputFolder,
                                         "gParametrized");
            volScalarField& T = exampleSteady._T();
            ITHACAstream::exportSolution(T,
                                         std::to_string(realTime),
                                         outputFolder,
                                         "T");
        }
        ITHACAstream::exportMatrix(timeVector, "time", "eigen", outputFolder);
        exampleSteady.postProcess(outputFolder, example.samplingTime, "gParametrized",
                                  example.probe1, example.probe2);
    }

    if (parameterizedBC_unsteadyTest)
    {
        word outputFolder = "./ITHACAoutput/testInverse/";
        example.assignTrueIF();
        example.set_gParametrized("rbf", shapeParameter, "linear");
        example.parameterizedBCoffline();
        volScalarField initialField = example.Ttrue[0];
        //volScalarField initialField = example._T();
        //scalar T0 = 650;
        //ITHACAutilities::assignIF(initialField, T0);
        example.parameterizedBC(outputFolder, initialField, linSysSolver, TSVDtruncation);
        example.inverseProblemPostProcess(outputFolder);
    }

    return 0;
}

