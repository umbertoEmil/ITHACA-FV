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
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "reducedInverseLaplacian.H"
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

#include "DanieliMold.H"
using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 0; //No verbose output
    double time;
    DanieliMold example(argc, argv);

    // Reading tests to perform
    ITHACAparameters* para= ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label TSVDtest = para->ITHACAdict->lookupOrDefault<int>("TSVDtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCerrorTest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest", 0);
    bool experimentalData =
        para->ITHACAdict->lookupOrDefault<bool>("experimentalData", 0);
    bool userDefinedTemperatures =
        para->ITHACAdict->lookupOrDefault<bool>("userDefinedTemperatures", 1);
    double moldWidth = 
        para->ITHACAdict->lookupOrDefault<double>("moldWidth", 1.97);
    double castWidth =
        para->ITHACAdict->lookupOrDefault<double>("castWidth", 0);

    // Reading para->eters from ITHACAdict
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


    double Tf = para->ITHACAdict->lookupOrDefault<double>("Tf", 300.0);
    double refGrad = para->ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction = para->ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);
    
    if(userDefinedTemperatures)
    {
        example.castWidth = moldWidth;
        // Performing true solution
        auto t1 = std::chrono::high_resolution_clock::now();
        example.solveTrue();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6;
        std::cout << "True solution took  = " << duration << " seconds" << std::endl;
        
        // Setting up the thermocouples
        example.readThermocouples();
        //example.getThermocouplesPlane();
        //Info << "Exiting" << endl;
        //exit(10);
        example.Tmeas.resize(example.thermocouplesCellID.size());
        example.Tdirect = example.Tmeas;
        example.Tsens = example.Tmeas;
        volScalarField& T(example._T());
        forAll(example.thermocouplesCellID, cellI)
        {
            example.Tmeas(cellI) = T.internalField()[example.thermocouplesCellID[cellI]];
        }
    }
    else if(experimentalData)
    {
        M_Assert(castWidth > 1e-6, "Set the castWidth");
        example.castWidth = castWidth;
        example.setBC();
        word TCfile =
            para->ITHACAdict->lookupOrDefault<word>("thermocouplesFile", "error");
        Info << "Reading thermocouples temperatures from file: " << TCfile << endl;
	M_Assert( TCfile.std::string::compare("error") != 1, "No thermocouples measurements file defined");
        example.readThermocouples();
	example.readTCfromFile(TCfile);
    }
    //
    //if (example.interpolation)
    //{
    //    Info << "Interpolating thermocouples measurements in the " <<
    //         "plane defined by the thermocouples" << endl;
    //    example.thermocouplesInterpolation();
    //}
    //else
    //{
    //    Info << "NOT interpolating thermocouples measurements" << endl;
    //}
    
    
    // Solving the inverse problem
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if (CGtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the CG inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/CGtest/";
        if(userDefinedTemperatures)
        {
            volScalarField gTrueField = example.list2Field(example.gTrue);
            ITHACAstream::exportSolution(gTrueField,
                                         "1", outputFolder,
                                         "gTrue");
	}

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
    
    if(parameterizedBCtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the parameterized BC inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCtest/";
        if(userDefinedTemperatures)
        {
            volScalarField gTrueField = example.list2Field(example.gTrue);
            ITHACAstream::exportSolution(gTrueField,
                                         "1", outputFolder,
                                         "gTrue");
	}

        List<word> linSys_solvers;
        linSys_solvers.resize(5);
        linSys_solvers[0] = "fullPivLU";
        linSys_solvers[1] = "jacobiSvd";
        linSys_solvers[2] = "householderQr";
        linSys_solvers[3] = "ldlt";
        linSys_solvers[4] = "TSVD";
        Eigen::VectorXd residualNorms;
        residualNorms.resize(linSys_solvers.size());
        example.set_gParametrized("rbf",0.9);
        example.parameterizedBCoffline();
    
        forAll(linSys_solvers, solverI)
        {
            Info << "Solver " << linSys_solvers[solverI] << endl;
            Info << endl;
    
            auto t1 = std::chrono::high_resolution_clock::now();
            example.parameterizedBC(linSys_solvers[solverI], 14);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6;
            std::cout << "Duration online part = " << duration << " seconds" << std::endl;
            volScalarField gParametrizedField = example.list2Field(example.g);
	    volScalarField& T(example._T());
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "gParametrized");
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
        if(userDefinedTemperatures)
        {
            example.postProcess(outputFolder, "gParametrized", userDefinedTemperatures);
	}
	else
	{
	    example.postProcess(outputFolder, "gParametrized", userDefinedTemperatures);
	}
        Info << "*********************************************************" << endl;
        Info << endl;
    
    }
    
    if(TSVDtest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the parameterized BC inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/TSVDtest/";
        if(userDefinedTemperatures)
        {
            volScalarField gTrueField = example.list2Field(example.gTrue);
            ITHACAstream::exportSolution(gTrueField,
                                         "1", outputFolder,
                                         "gTrue");
	}

        Eigen::VectorXd residualNorms;
        example.set_gParametrized("rbf",0.9);
        example.parameterizedBCoffline();
    
        Eigen::VectorXd tsvdTruncation = Eigen::VectorXd::LinSpaced(20, 10, 30);
        residualNorms.resize(tsvdTruncation.size());
        for(int i = 0; i < tsvdTruncation.size(); i++)
        {
    
            example.parameterizedBC("TSVD", tsvdTruncation(i));
	    
            volScalarField gParametrizedField = example.list2Field(example.g);
	    volScalarField& T(example._T());
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "gParametrized");
            ITHACAstream::exportSolution(T,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "T");
            residualNorms(i) = Foam::sqrt(
                        example.residual.squaredNorm());
        }
        Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
    
        ITHACAstream::exportVector(residualNorms, "residuals2norm", "eigen",
                                   outputFolder);
        example.postProcess(outputFolder, "gParametrized", userDefinedTemperatures);
        Info << "*********************************************************" << endl;
        Info << endl;
    
    }
    //
    //    if(CGnoiseTest)
    //    {
    //        Info << endl;
    //        Info << "*********************************************************" << endl;
    //        Info << "Testing CG inverse solver with NOISY data" << endl;
    //        Info << "Performing " << Ntests << " tests." << endl;
    //        Info << endl;
    //        word outputFolder = "./ITHACAoutput/CGnoiseTest/";
    //  example.gTrue = -example.heatFlux_hotSide / example.k;
    //  volScalarField gTrueField = example.list2Field(example.gTrue);
    //        ITHACAstream::exportSolution(gTrueField,
    //                                     "1", outputFolder,
    //                                     "gTrue");
    //        example.saveSolInLists = 1;
    //  Eigen::VectorXd variance = example.Tmeas * 0.02;
    //  std::cout  << "variance = " << variance << std::endl;
    //  std::cout  << "example.Tmeas = " << example.Tmeas << std::endl;
    //
    //  Info << "Noise variance mean = " << variance.mean() << endl;
    //        Eigen::VectorXd TmeasOrig = example.Tmeas;
    //  Info << "I use the discrepancy principle as stopping criterium for CG" << endl;
    //  example.Jtol = example.Tmeas.size() * variance.maxCoeff() * variance.maxCoeff() / 2;
    //  Info << "Stopping for J < " << example.Jtol << endl;
    //
    //        for (label i = 0; i < Ntests; i++)
    //        {
    //            Info << "Test " << i << endl;
    //            example.Tmeas = TmeasOrig;
    //            Eigen::VectorXd measurementsError(example.Tmeas.size());
    //            for(int i = 0; i < example.Tmeas.size(); i++)
    //            {
    //                measurementsError(i) = variance(i) * stochastic::set_normal_random(0.0, 1.0);
    //            }
    //            example.Tmeas += measurementsError;
    //            Info << "Measurements error L2 norm= " << measurementsError.norm() <<
    //                 endl;
    //            Info << "example.gList.size() = " << example.gList.size() << endl;
    //      if (example.conjugateGradient())
    //            {
    //                Info << "CG converged" << endl;
    //                volScalarField heatFluxField = example.list2Field(example.gList[example.gList.size() - 1]);
    //                ITHACAstream::exportSolution(heatFluxField,
    //                                             std::to_string(i + 1), outputFolder,
    //                                             "g_CG");
    //                Info << "************************************" << endl;
    //      Info << endl << endl;
    //            }
    //            else
    //            {
    //                Info << "CG did not converged" << endl;
    //                Info << "************************************" << endl;
    //      Info << endl << endl;
    //            }
    //
    //        }
    //
    //        example.postProcess(outputFolder, "g_CG");
    //
    //        Info << "*********************************************************" << endl;
    //        Info << endl;
    //    }
    //
    //    if (para->eterizedBCerrorTest)
    //    {
    //        Info << endl;
    //        Info << "*********************************************************" << endl;
    //        Info << "Testing para->eterized BC inverse solver with NOISY data" <<
    //             endl;
    //        word outputFolder = "./ITHACAoutput/para->eterizedBCerrorTest/";
    //  example.gTrue = -example.heatFlux_hotSide / example.k;
    //  volScalarField gTrueField = example.list2Field(example.gTrue);
    //        ITHACAstream::exportSolution(gTrueField,
    //                                     "1", outputFolder,
    //                                     "gTrue");
    //        List<List<scalar>> heatFluxWeights;
    //        Eigen::VectorXd residualNorms;
    //        scalar innerField = 1.0;
    //        example.set_gParametrized("rbf", 0.7);
    //        example.para->eterizedBCoffline();
    //        List<word> linSys_solvers;
    //        linSys_solvers.resize(1);
    //        linSys_solvers[0] = "TSVD";
    //        //linSys_solvers[0] = "fullPivLU";
    //        //linSys_solvers[1] = "jacobiSvd";
    //        //linSys_solvers[2] = "householderQr";
    //        //linSys_solvers[3] = "ldlt";
    //        //linSys_solvers[4] = "TSVD";
    //        Info << "Introducing error in the measurements" << endl;
    //        Info << "Performing " << Ntests << " tests." << endl;
    //        residualNorms.resize(Ntests * linSys_solvers.size());
    //        Eigen::VectorXd TmeasOrig = example.Tmeas;
    //
    //        for (label i = 0; i < Ntests; i++)
    //        {
    //            Info << "Test " << i << endl;
    //            example.Tmeas = TmeasOrig;
    //            Eigen::VectorXd measurementsError(example.Tmeas.size());
    //            for(int i = 0; i < example.Tmeas.size(); i++)
    //            {
    //                measurementsError(i) = example.Tmeas.mean() * 0.02 * stochastic::set_normal_random(0.0, 1.0);
    //            }
    //            example.Tmeas += measurementsError;
    //
    //            List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
    //            List<scalar> solutionNorms;
    //            forAll(linSys_solvers, solverI)
    //            {
    //                Info << "Solver " << linSys_solvers[solverI] << endl;
    //                example.para->eterizedBC(outputFolder, linSys_solvers[solverI], 3);
    //                Info << endl;
    //                volScalarField gParametrizedField = example.list2Field(example.g);
    //                ITHACAstream::exportSolution(gParametrizedField,
    //                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
    //                                             outputFolder,
    //                                             "gParametrized");
    //                ITHACAstream::exportSolution(example.T,
    //                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
    //                                             outputFolder,
    //                                             "T");
    //                residualNorms(i * linSys_solvers.size() + solverI) = Foam::sqrt(
    //                            example.residual.squaredNorm());
    //            }
    //            Info << "Measurements error L2 norm= " << measurementsError.norm() <<
    //                 endl;
    //        }
    //
    //        example.postProcess(outputFolder,"gParametrized", innerField);
    //        Info << "*********************************************************" << endl;
    //        Info << endl;
    //    }
    return 0;
}

