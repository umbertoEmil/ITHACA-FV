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
    totalHeatFluxTest.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem_CG.H"
#include "inverseLaplacianProblem_paramBC.H"
#include "inverseLaplacianProblemTotalHeatMeasure_CG.H"
#include "inverseLaplacianProblemTotalHeatMeasure_paramBC.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "totalHeatFluxTest.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    double time;
    totalHeatFluxTest_CG example_CG(argc, argv);
    totalHeatFluxTest_paramBC example_paramBC(argc, argv);
    // Reading parameters from ITHACAdict
    ITHACAparameters* para = ITHACAparameters::getInstance(example_CG._mesh(),
                             example_CG._runTime());
    // Reading parameters from ITHACAdict
    example_CG.g_0 = para->ITHACAdict->lookupOrDefault<double>("g_0", 0.0); 
    example_CG.g_X = para->ITHACAdict->lookupOrDefault<double>("g_X", 0.0); 
    example_CG.g_Z = para->ITHACAdict->lookupOrDefault<double>("g_Z", 0.0); 
    example_CG.Tf_0 = para->ITHACAdict->lookupOrDefault<double>("Tf_0", 0.0); 
    example_CG.Tf_delta = para->ITHACAdict->lookupOrDefault<double>("Tf_delta", 0.0); 

    example_paramBC.g_0 = example_CG.g_0;
    example_paramBC.g_X = example_CG.g_X;
    example_paramBC.g_Z = example_CG.g_Z;
    example_paramBC.Tf_0 = example_CG.Tf_0;
    example_paramBC.Tf_delta = example_CG.Tf_delta;

    word solver = para->ITHACAdict->lookupOrDefault<word>("linSysSolver", "fullPivLU");
    int TSVDtrunc = para->ITHACAdict->lookupOrDefault<int>("TSVDregularization", 0);

    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);

    // Reading parameters from ITHACAdict
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);

    example_CG.interpolation = para->ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                    0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);
    
    example_CG.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_CG.k > 0, "thermalConductivity, k, not specified");
    example_CG.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_CG.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.gIntegralWeight = para->ITHACAdict->lookupOrDefault<double>("gIntegralWeight", 0);
    example_paramBC.gIntegralWeight = example_CG.gIntegralWeight;
    example_paramBC.k = example_CG.k;
    example_paramBC.H = example_CG.H;

    //*******************************************//
    fvMesh& mesh = example_CG._mesh();
    example_CG.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = mesh.boundaryMesh()[example_CG.hotSide_ind].size(); 
    example_CG.g.resize(hotSideSize);
    example_CG.gTrue.resize(hotSideSize);
    
    forAll(example_CG.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example_CG.hotSide_ind].faceCentres()[faceI].x();
        scalar faceZ =
            mesh.boundaryMesh()[example_CG.hotSide_ind].faceCentres()[faceI].z();
        example_CG.g[faceI] = (example_CG.g_X * (faceX - 1) * (faceX - 1) + example_CG.g_Z * faceZ + example_CG.g_0) ;
    }

    example_CG.set_Tf();
    example_CG.gTrue = example_CG.g;
    example_paramBC.gTrue = example_CG.g;
    example_CG.solveTrue();
    volScalarField& T(example_CG._T());
    example_paramBC.set_Tf();
    example_paramBC.g = example_CG.g;
    example_paramBC.solveTrue();


    // Setting up the thermocouples
    example_CG.gIntegral_meas = ITHACAutilities::integralOnPatch(mesh, example_CG.gTrue, "hotSide");
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T);
    std::cout << "debug: Tmeas = " << example_CG.Tmeas << std::endl;
    example_paramBC.gIntegral_meas = example_CG.gIntegral_meas;
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = example_CG.Tmeas;

    // Solving the inverse problem
    if (CGtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the CG inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/CGtest/";
        volScalarField gTrueField = example_CG.list2Field(example_CG.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        example_CG.saveSolInLists = 1;
        auto t1 = std::chrono::high_resolution_clock::now();

        if (example_CG.conjugateGradient())
        {
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                            ( t2 - t1 ).count() / 1e6;
            std::cout << "Duration = " << duration << " seconds" << std::endl;
            Info << "CG converged" << endl;
            PtrList<volScalarField> heatFluxField;
            forAll(example_CG.gList, solutionI)
            {
                heatFluxField.append(example_CG.list2Field(example_CG.gList[solutionI]).clone());
                ITHACAstream::exportSolution(heatFluxField[solutionI],
                                             std::to_string(solutionI + 1), outputFolder,
                                             "g_CG");
            }
            example_CG.postProcess(outputFolder, "g_CG");
        }
        else
        {
            Info << "CG did not converged" << endl;
        }

        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (parameterizedBCtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the parameterized BC inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCtest/";
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<word> linSys_solvers;
        linSys_solvers.resize(5);
        linSys_solvers[0] = "fullPivLU";
        linSys_solvers[1] = "jacobiSvd";
        linSys_solvers[2] = "householderQr";
        linSys_solvers[3] = "ldlt";
        linSys_solvers[4] = "TSVD";
        Eigen::VectorXd residualNorms;
        residualNorms.resize(linSys_solvers.size());
        example_paramBC.set_gParametrized("rbf", 0.7);
        example_paramBC.parameterizedBCoffline();
        forAll(linSys_solvers, solverI)
        {
            Info << "Solver " << linSys_solvers[solverI] << endl;
            Info << endl;
            auto t1 = std::chrono::high_resolution_clock::now();
            example_paramBC.parameterizedBC(linSys_solvers[solverI], TSVDtrunc);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                            ( t2 - t1 ).count() / 1e6;
            std::cout << "Duration online part = " << duration << " seconds" << std::endl;
            volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "gParametrized");
	    volScalarField& T = example_paramBC._T();
            ITHACAstream::exportSolution(T,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "T");
            residualNorms(solverI) = Foam::sqrt(
                                         example_paramBC.residual.squaredNorm());
        }
        Eigen::MatrixXd A = example_paramBC.Theta.transpose() * example_paramBC.Theta + example_paramBC.gIntegralWeight * example_paramBC.Phi;
        ITHACAstream::exportMatrix(residualNorms, "residuals2norm", "eigen",
                                   outputFolder);
        example_paramBC.postProcess(outputFolder, "gParametrized");
        Info << "*********************************************************" << endl;
        Info << endl;
    }
    return 0;
}


