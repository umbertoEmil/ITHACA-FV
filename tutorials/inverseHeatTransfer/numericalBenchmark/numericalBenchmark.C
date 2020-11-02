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
    numericalBenchmark.C
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
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "numericalBenchmark.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"


using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    double time;
    numericalBenchmark_CG example_CG(argc, argv);
    numericalBenchmark_paramBC example_paramBC(argc, argv);
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
    label TSVDtruncation = para->ITHACAdict->lookupOrDefault<label>("TSVDtruncation", 3);

    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCtest_RBFwidth =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest_RBFwidth", 0);
    label parameterizedBCerrorTest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest", 0);
    label parameterizedBCerrorTest_TSVD = 
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest_TSVD", 0);



    // Reading parameters from ITHACAdict
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example_CG.thermocouplesNum =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    example_paramBC.thermocouplesNum = example_CG.thermocouplesNum;
    M_Assert(example_CG.thermocouplesNum > 0, "Number of thermocouples not specified");

    example_CG.interpolation = para->ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                    0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);

    int TSVDtrunc = para->ITHACAdict->lookupOrDefault<int>("TSVDregularization", 0);

    double noiseLevel = para->ITHACAdict->lookupOrDefault<double>("noiseLevel",
                   0);
    example_CG.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_CG.k > 0, "thermalConductivity, k, not specified");
    example_CG.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    example_paramBC.k = example_CG.k;
    example_paramBC.H = example_CG.H;
    M_Assert(example_CG.H > 0, "Heat transfer coeff, H, not specified");
    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double refGrad = para->ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction = para->ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);

    //*******************************************//
    fvMesh& mesh = example_CG._mesh();
    example_CG.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = mesh.boundaryMesh()[example_CG.hotSide_ind].size(); 
    example_CG.g.resize(hotSideSize);
    example_CG.gTrue.resize(hotSideSize);
    double x0 = 0.788;
    double z0 = 0.6;
    double radius = 0.2;
    
    forAll(example_CG.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example_CG.hotSide_ind].faceCentres()[faceI].x();
        scalar faceZ =
            mesh.boundaryMesh()[example_CG.hotSide_ind].faceCentres()[faceI].z();
	//example_CG.g[faceI] = - example_CG.k * (example_CG.b * std::exp(- 5 * ((faceX - x0) * (faceX - x0) + (faceZ - z0) * (faceZ - z0))));
	//if( (faceX - x0) * (faceX - x0) + (faceZ - z0) * (faceZ - z0) < radius * radius)
	//{
	//   example_CG.g[faceI] = - example_CG.k * (example_CG.b);
	//}
	//else 
	//{
	//   example_CG.g[faceI] = 0;
	//}
        example_CG.g[faceI] = (example_CG.g_X * (faceX - 1) * (faceX - 1) + example_CG.g_Z * faceZ + example_CG.g_0) ;
    }

    example_CG.set_Tf();
    example_CG.gTrue = example_CG.g;
    example_CG.solveTrue();

    Info << "debug 1" << endl;
    Info << "debug 2" << endl;
    example_paramBC.g = example_CG.g;
    example_paramBC.gTrue = example_CG.gTrue;
    example_paramBC.set_Tf();
    example_paramBC.solveTrue();
    Info << "debug 3" << endl;
    
    volScalarField& T(example_CG._T());

    // Setting up the thermocouples
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T);
    std::cout << "debug: Tmeas = " << example_CG.Tmeas << std::endl;
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = example_paramBC.fieldValueAtThermocouples(T);


    if (example_CG.interpolation)
    {
        Info << "Interpolating thermocouples measurements in the " <<
             "plane defined by the thermocouples" << endl;
        example_CG.thermocouplesInterpolation();
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
                heatFluxField.append(example_CG.list2Field(example_CG.gList[solutionI]));
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
            example_paramBC.parameterizedBC(linSys_solvers[solverI], TSVDtrunc);
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
        Eigen::MatrixXd A = example_paramBC.Theta.transpose() * example_paramBC.Theta;
        ITHACAstream::exportMatrix(residualNorms, "residuals2norm", "eigen",
                                   outputFolder);
        example_paramBC.postProcess(outputFolder, "gParametrized");
        Info << "*********************************************************" << endl;
        Info << endl;
    }



    if (parameterizedBCtest_RBFwidth)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Performing test for the parameterized BC inverse solver" << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCtest_RBFparameter/";
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        
        int rbfWidth_size = 9;
        Eigen::VectorXd rbfWidth(rbfWidth_size);
        rbfWidth << 100, 10, 1, 0.33, 0.1, .033, 0.01, 0.0033, 0.001;

        Eigen::VectorXd residualNorms;
        residualNorms.resize(rbfWidth_size);
        Eigen::VectorXd heatFluxL2norm(rbfWidth_size);
        Eigen::VectorXd heatFluxLinfNorm = heatFluxL2norm;
        Eigen::VectorXd condNumber = heatFluxL2norm;
        Eigen::MatrixXd singVal;

        for(int i = 0; i < rbfWidth_size; i++)
        {
            Info << "*********************************************************" << endl;
            Info << "RBF parameter " << rbfWidth(i) << endl;
            Info << endl;

            example_paramBC.set_gParametrized("rbf", rbfWidth(i));
            example_paramBC.parameterizedBCoffline(1);
            example_paramBC.parameterizedBC("fullPivLU");

            volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "gParametrized");
	    volScalarField& T(example_paramBC._T());
            ITHACAstream::exportSolution(T,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "T");
            residualNorms(i) = Foam::sqrt(
                                         example_paramBC.residual.squaredNorm());
            Eigen::MatrixXd A = example_paramBC.Theta.transpose() * example_paramBC.Theta;
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::MatrixXd singularValues = svd.singularValues();
            singVal.conservativeResize(singularValues.rows(), singVal.cols() + 1);
            singVal.col(i) = singularValues;
            double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
            Info << "Condition number = " << conditionNumber << endl;
            condNumber(i) = conditionNumber;


            volScalarField gDiffField = gParametrizedField - gTrueField;
            scalar EPS = 1e-6;
            volScalarField relativeErrorField(gTrueField);
            for (label i = 0; i < relativeErrorField.internalField().size(); i++)
            {
                if (std::abs(gTrueField.ref()[i]) < EPS)
                {
                    relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
                }
                else
                {
                    relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / gTrueField.ref()[i];
                }
            }
            ITHACAstream::exportSolution(relativeErrorField,
                                         std::to_string(i + 1), outputFolder,
                                         "relativeErrorField");
            heatFluxL2norm(i) = ITHACAutilities::L2normOnPatch(mesh, relativeErrorField,
                                        "hotSide");
            heatFluxLinfNorm(i) = ITHACAutilities::LinfNormOnPatch(mesh, relativeErrorField,
                                        "hotSide");
        }

        ITHACAstream::exportMatrix(condNumber, "condNumber", "eigen",
                                   outputFolder);
        ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                                   outputFolder);
        ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                                   outputFolder);
        ITHACAstream::exportMatrix(singVal, "singularValues", "eigen",
                                   outputFolder);
        ITHACAstream::exportMatrix(residualNorms, "residuals2norm", "eigen",
                                   outputFolder);
        example_paramBC.postProcess(outputFolder, "gParametrized");
        Info << "*********************************************************" << endl;
        Info << "*********************************************************" << endl;
        Info << endl;
    }




    if (parameterizedBCerrorTest_TSVD)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing parameterized BC TSVD with NOISY data" <<
             endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCnoiseTest_TSVD/";
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<List<scalar>> heatFluxWeights;
        scalar innerField = 1.0;
        example_paramBC.set_gParametrized("rbf", 0.1);
        example_paramBC.parameterizedBCoffline();

        Info << "Introducing error in the measurements" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;

        Eigen::VectorXi TSVDtruc(15);
        TSVDtruc << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14, 15;
        std::cout << "debug : TSVDtruc = \n" << TSVDtruc.transpose() << std::endl;
        std::cout << "debug : Tmeas orig = \n" << example_paramBC.Tmeas.transpose() << std::endl;

        int Ntrunc = TSVDtruc.size();

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example_paramBC.addNoise(noiseLevel);
            std::cout << "debug : Tmeas = \n" << example_paramBC.Tmeas.transpose() << std::endl;

            List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
            List<scalar> solutionNorms;
            for(int truncI = 0; truncI < Ntrunc; truncI++)
            {
                example_paramBC.parameterizedBC("TSVD", TSVDtruc(truncI));
                volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(i * Ntrunc + truncI + 1),
                                             outputFolder,
                                             "gParametrized");
                ITHACAstream::exportSolution(example_paramBC.T,
                                             std::to_string(i * Ntrunc + truncI + 1),
                                             outputFolder,
                                             "T");
            }
            example_paramBC.Tmeas = TmeasOrig;
        }

        example_paramBC.postProcess(outputFolder, "gParametrized", innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    return 0;
}


