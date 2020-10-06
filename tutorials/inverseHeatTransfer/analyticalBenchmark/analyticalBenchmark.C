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
    analyticalBenchmark.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "interpolation.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "reducedInverseLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "analyticalBenchmark.H"


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    analyticalBenchmark example(argc, argv);
    //Setting parameters for the analytical benchmark
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    example.a = a;
    example.b = b;
    example.c = c;
    example.d = d;

    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCerrorTest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest", 0);
    // Reading parameters from ITHACAdict
    example.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example.thermocouplesNum =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    M_Assert(example.thermocouplesNum > 0, "Number of thermocouples not specified");
    example.interpolation = para->ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                    0.000001);
    example.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);
    example.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);

    // setting analytical solution
    volScalarField T_true(example._T());

    for (label i = 0; i < T_true.internalField().size(); i++)
    {
        auto cx = T_true.mesh().C()[i].component(vector::X);
        auto cy = T_true.mesh().C()[i].component(vector::Y);
        auto cz = T_true.mesh().C()[i].component(vector::Z);
        T_true.ref()[i] = a * cx * cx + b * cx * cy + c * cy - a * cz * cz + c;
    }
    //da correggere, si é aggiunto il termine c * cy 
    //scalar Ttrue_L2norm = std::sqrt(28 / 45 * a * a + 1 / 9 * b * b + c * c + 3 * a * b / 8 + 4/3 *a * c + b * c / 2) ;
    //scalar Ttrue_L2norm = std::sqrt(13 / 45 * a * a + 1 / 9 * b * b + 7/3 * c * c + 1 / 12 * a * b - 1 / 3 *a * c + 5 / 12 * b * c) ;
    scalar Ttrue_L2norm = 25.8789694352;


    fvMesh& mesh = example._mesh();
    example.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T_true.boundaryField()[example.hotSide_ind].size();
    example.g.resize(hotSideSize);
    example.gTrue.resize(hotSideSize);
    forAll(example.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example.hotSide_ind].faceCentres()[faceI].x();
        example.g[faceI] = example.k * (b * faceX + c) ;
    }
    example.gTrue = example.g; //-example.heatFlux_hotSide / example.k;
    example.solveTrue();
    

    volScalarField& T(example._T());
    Info << "Exporting analytical solution" << endl;
    ITHACAstream::exportSolution(T_true, "1", "./ITHACAoutput/true/",
                                 "analyticalSol");
    volScalarField error = T_true - T;
    ITHACAstream::exportSolution(error, "1", "./ITHACAoutput/true/", "error");
    Info << "L2 norm of the relative error = " << ITHACAutilities::errorL2Rel(
             T_true, T) << endl;
    Info << "Linf norm of the relative error = " <<
         ITHACAutilities::errorLinfRel(T_true, T) << endl;
    Info << "L2 norm of the error = " << ITHACAutilities::errorL2Abs(T_true, T) << endl;
    Info << "Linf norm of the error = " << ITHACAutilities::LinfNorm(error) << endl;
    Info << "L2 norm of T = " << ITHACAutilities::L2Norm(T) << endl;
    Info << "L2 norm of Tanal = " << Ttrue_L2norm << endl;
    Info << "L2 norm of diff = " << ITHACAutilities::L2Norm(T) - Ttrue_L2norm  <<  endl;
    Info << "L2 norm of rel diff = " << (ITHACAutilities::L2Norm(T) - Ttrue_L2norm) / Ttrue_L2norm  << endl<< endl;
    // Setting up the thermocouples
    example.readThermocouples();
    example.Tmeas = example.fieldValueAtThermocouples(T_true);
    std::cout << "debug: Tmeas = " << example.Tmeas << std::endl;

    // Introducing error in the measurements
    //Tmeas += ITHACAutilities::rand(Tmeas.size(), 1, -2, 2);
    //Eigen::VectorXd measurementsError(Tmeas.size());
    //for(int i = 0; i < Tmeas.size(); i++)
    //{
    //    measurementsError(i) = Tmeas.mean() * 0.02 * stochastic::set_normal_random(0.0, 1.0);
    //}
    //Tmeas += measurementsError;

    if (example.interpolation)
    {
        Info << "Interpolating thermocouples measurements in the " <<
             "plane defined by the thermocouples" << endl;
        example.thermocouplesInterpolation();
    }
    else
    {
        Info << "NOT interpolating thermocouples measurements" << endl;
    }

    // Solving the inverse problem
    if (CGtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the CG inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/CGtest/";
        volScalarField gTrueField = example.list2Field(example.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        example.saveSolInLists = 1;
        auto t1 = std::chrono::high_resolution_clock::now();

        if (example.conjugateGradient())
        {
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                            ( t2 - t1 ).count() / 1e6;
            std::cout << "Duration = " << duration << " seconds" << std::endl;
            Info << "CG converged" << endl;
            PtrList<volScalarField> heatFluxField;
            forAll(example.gList, solutionI)
            {
                heatFluxField.append(example.list2Field(example.gList[solutionI]));
                ITHACAstream::exportSolution(heatFluxField[solutionI],
                                             std::to_string(solutionI + 1), outputFolder,
                                             "g_CG");
            }
            example.postProcess(outputFolder, "g_CG");
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
        volScalarField gTrueField = example.list2Field(example.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<word> linSys_solvers;
        linSys_solvers.resize(1);
        linSys_solvers[0] = "fullPivLU";
        //linSys_solvers[1] = "jacobiSvd";
        //linSys_solvers[2] = "householderQr";
        //linSys_solvers[3] = "ldlt";
        //linSys_solvers[4] = "TSVD";
        Eigen::VectorXd residualNorms;
        residualNorms.resize(linSys_solvers.size());
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        forAll(linSys_solvers, solverI)
        {
            Info << "Solver " << linSys_solvers[solverI] << endl;
            Info << endl;
            auto t1 = std::chrono::high_resolution_clock::now();
            example.parameterizedBC(linSys_solvers[solverI], 3);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                            ( t2 - t1 ).count() / 1e6;
            std::cout << "Duration online part = " << duration << " seconds" << std::endl;
            volScalarField gParametrizedField = example.list2Field(example.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "gParametrized");
	    volScalarField& T(example._T());
            ITHACAstream::exportSolution(T,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "T");
            residualNorms(solverI) = Foam::sqrt(
                                         example.residual.squaredNorm());
        }
        Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
        ITHACAstream::exportVector(residualNorms, "residuals2norm", "eigen",
                                   outputFolder);
        example.postProcess(outputFolder, "gParametrized");
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (CGnoiseTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing CG inverse solver with NOISY data" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/CGnoiseTest/";
        example.gTrue = -example.heatFlux_hotSide / example.k;
        volScalarField gTrueField = example.list2Field(example.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        example.saveSolInLists = 1;
        Eigen::VectorXd variance = example.Tmeas * 0.02;
        std::cout  << "variance = " << variance << std::endl;
        std::cout  << "example.Tmeas = " << example.Tmeas << std::endl;
        Info << "Noise variance mean = " << variance.mean() << endl;
        Eigen::VectorXd TmeasOrig = example.Tmeas;
        Info << "I use the discrepancy principle as stopping criterium for CG" << endl;
        example.Jtol = example.Tmeas.size() * variance.maxCoeff() *
                       variance.maxCoeff() / 2;
        Info << "Stopping for J < " << example.Jtol << endl;

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example.Tmeas = TmeasOrig;
            Eigen::VectorXd measurementsError(example.Tmeas.size());

            for (int i = 0; i < example.Tmeas.size(); i++)
            {
                //measurementsError(i) = variance(i) * stochastic::set_normal_random(0.0, 1.0);
		Info << "add stochastic.H" << endl;
		exit(10);
			
            }

            example.Tmeas += measurementsError;
            Info << "Measurements error L2 norm= " << measurementsError.norm() <<
                 endl;
            Info << "example.gList.size() = " << example.gList.size() << endl;

            if (example.conjugateGradient())
            {
                Info << "CG converged" << endl;
                volScalarField heatFluxField = example.list2Field(
                                                   example.gList[example.gList.size() - 1]);
                ITHACAstream::exportSolution(heatFluxField,
                                             std::to_string(i + 1), outputFolder,
                                             "g_CG");
                Info << "************************************" << endl;
                Info << endl << endl;
            }
            else
            {
                Info << "CG did not converged" << endl;
                Info << "************************************" << endl;
                Info << endl << endl;
            }
        }

        example.postProcess(outputFolder, "g_CG");
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (parameterizedBCerrorTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing parameterized BC inverse solver with NOISY data" <<
             endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCerrorTest/";
        example.gTrue = -example.heatFlux_hotSide / example.k;
        volScalarField gTrueField = example.list2Field(example.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<List<scalar>> heatFluxWeights;
        Eigen::VectorXd residualNorms;
        scalar innerField = 1.0;
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        List<word> linSys_solvers;
        linSys_solvers.resize(1);
        linSys_solvers[0] = "TSVD";
        //linSys_solvers[0] = "fullPivLU";
        //linSys_solvers[1] = "jacobiSvd";
        //linSys_solvers[2] = "householderQr";
        //linSys_solvers[3] = "ldlt";
        //linSys_solvers[4] = "TSVD";
        Info << "Introducing error in the measurements" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        residualNorms.resize(Ntests * linSys_solvers.size());
        Eigen::VectorXd TmeasOrig = example.Tmeas;

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example.Tmeas = TmeasOrig;
            Eigen::VectorXd measurementsError(example.Tmeas.size());

            for (int i = 0; i < example.Tmeas.size(); i++)
            {
                //measurementsError(i) = example.Tmeas.mean() * 0.02 *
                //                       stochastic::set_normal_random(0.0, 1.0);
		Info << "add stochastic.H" << endl;
		exit(10);
            }

            example.Tmeas += measurementsError;
            List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
            List<scalar> solutionNorms;
            forAll(linSys_solvers, solverI)
            {
                Info << "Solver " << linSys_solvers[solverI] << endl;
                example.parameterizedBC(linSys_solvers[solverI], 3);
                Info << endl;
                volScalarField gParametrizedField = example.list2Field(example.g);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
                                             outputFolder,
                                             "gParametrized");
                ITHACAstream::exportSolution(example.T,
                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
                                             outputFolder,
                                             "T");
                residualNorms(i * linSys_solvers.size() + solverI) = Foam::sqrt(
                            example.residual.squaredNorm());
            }
            Info << "Measurements error L2 norm= " << measurementsError.norm() <<
                 endl;
        }

        example.postProcess(outputFolder, "gParametrized", innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    return 0;
}

