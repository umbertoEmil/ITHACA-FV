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
    example_paramBC.of a heat transfer Reduction Problem
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
#include "inverseLaplacianProblem_CG.H"
#include "inverseLaplacianProblem_paramBC.H"
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
    analyticalBenchmark_paramBC example_paramBC(argc, argv);
    analyticalBenchmark_CG example_CG(argc, argv);
    //Setting parameters for the analytical benchmark
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    example_paramBC.a = a;
    example_paramBC.b = b;
    example_paramBC.c = c;
    example_paramBC.d = d;

    example_CG.a = a;
    example_CG.b = b;
    example_CG.c = c;
    example_CG.d = d;

    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example_paramBC._mesh(),
                             example_paramBC._runTime());
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label CGnoiseLevelTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseLevelTest", 0);
    label ParamBCnoiseLevelTest = para->ITHACAdict->lookupOrDefault<int>("ParamBCnoiseLevelTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCtest_RBFwidth =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest_RBFwidth", 0);
    label parameterizedBCerrorTest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest", 0);
    label CGandParamErrorTest =
        para->ITHACAdict->lookupOrDefault<int>("CGandParamErrorTest", 0);
    label parameterizedBCerrorTest_TSVD =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest_TSVD", 0);
    // Reading parameters from ITHACAdict
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example_CG.interpolation = para->ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                    0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);

    example_paramBC.thermocouplesNum =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    example_CG.thermocouplesNum = example_paramBC.thermocouplesNum;
    M_Assert(example_paramBC.thermocouplesNum > 0, "Number of thermocouples not specified");


    example_paramBC.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_paramBC.k > 0, "thermalConductivity, k, not specified");
    example_paramBC.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_paramBC.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.k = example_paramBC.k;
    example_CG.H = example_paramBC.H;


    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double noiseLevel = para->ITHACAdict->lookupOrDefault<double>("noiseLevel",
                   0);

    // setting analytical solution
    volScalarField T_true(example_paramBC._T());

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


    fvMesh& mesh = example_paramBC._mesh();
    example_paramBC.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T_true.boundaryField()[example_paramBC.hotSide_ind].size();
    example_paramBC.g.resize(hotSideSize);
    example_paramBC.gTrue.resize(hotSideSize);
    forAll(example_paramBC.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example_paramBC.hotSide_ind].faceCentres()[faceI].x();
        example_paramBC.g[faceI] = example_paramBC.k * (b * faceX + c) ;
    }
    example_paramBC.gTrue = example_paramBC.g;
    example_paramBC.solveTrue();

    example_CG.g = example_paramBC.g;
    example_CG.gTrue = example_paramBC.g;
    example_CG.solveTrue();
    

    volScalarField& T(example_paramBC._T());
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
    example_paramBC.readThermocouples();
    example_paramBC.Tmeas = example_paramBC.fieldValueAtThermocouples(T_true);

    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T_true);

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
        linSys_solvers.resize(1);
        linSys_solvers[0] = "fullPivLU";
        Eigen::VectorXd residualNorms;
        residualNorms.resize(linSys_solvers.size());
        auto tO1 = std::chrono::high_resolution_clock::now();
        example_paramBC.set_gParametrized("rbf", 0.1);
        example_paramBC.parameterizedBCoffline();
        auto tO2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                        ( tO2 - tO1 ).count() / 1e6;
        std::cout << "Duration offline part = " << duration << " seconds" << std::endl;
        forAll(linSys_solvers, solverI)
        {
            Info << "Solver " << linSys_solvers[solverI] << endl;
            Info << endl;
            auto t1 = std::chrono::high_resolution_clock::now();
            example_paramBC.parameterizedBC(linSys_solvers[solverI], 6);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>
                            ( t2 - t1 ).count() / 1e6;
            std::cout << "Duration online part = " << duration << " seconds" << std::endl;
            volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(solverI + 1),
                                         outputFolder,
                                         "gParametrized");
	    volScalarField& T(example_paramBC._T());
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
        rbfWidth << 100, 10, 1, 0.33, 0.1, .033, 0.01, 0.001, 0.0001;

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

    if (CGnoiseTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing CG inverse solver with NOISY data" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        Info << endl;
        word outputFolder = "./ITHACAoutput/CGnoiseTest/";
        volScalarField gTrueField = example_CG.list2Field(example_CG.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        example_CG.saveSolInLists = 1;
        Eigen::VectorXd TmeasOrig = example_CG.Tmeas;

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example_CG.addNoise(noiseLevel);
            Info << "Stopping for J < " << example_CG.Jtol << endl;


            if (example_CG.conjugateGradient())
            {
                Info << "CG converged" << endl;
                volScalarField heatFluxField = example_CG.list2Field(
                                                   example_CG.gList[example_CG.gList.size() - 1]);
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
            example_CG.Tmeas = TmeasOrig;
        }

        example_CG.postProcess(outputFolder, "g_CG");
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (CGnoiseLevelTest || ParamBCnoiseLevelTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing CG inverse solver with NOISY data" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        Info << endl;
        word outputFolder;
        if(CGnoiseLevelTest)
        {
            outputFolder = "./ITHACAoutput/CGnoiseLevelTest/";
        }
        else if(ParamBCnoiseLevelTest)
        {
            outputFolder = "./ITHACAoutput/ParamBCnoiseLevelTest/";
            example_paramBC.set_gParametrized("rbf", 0.1);
            example_paramBC.parameterizedBCoffline();
        }
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        example_CG.saveSolInLists = 1;
        Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;

        Eigen::VectorXd noiseLevelVec(10);
        noiseLevelVec << .005, .01, .02, .03, .04, .05, .06, .07, .08, .1;


        for (label NLi = 0; NLi < noiseLevelVec.size() ; NLi++)
        {
            noiseLevel = noiseLevelVec(NLi);
            for (label i = 0; i < Ntests; i++)
            {
                Info << "Test " << i << endl;
                example_paramBC.addNoise(noiseLevel);
                example_CG.Tmeas = example_paramBC.Tmeas;
                Info << "Stopping for J < " << example_CG.Jtol << endl;


                if(CGnoiseLevelTest)
                {
                    if (example_CG.conjugateGradient())
                    {
                        Info << "CG converged" << endl;
                        volScalarField heatFluxField = example_CG.list2Field(
                                                           example_CG.gList[example_CG.gList.size() - 1]);
                        ITHACAstream::exportSolution(heatFluxField,
                                                     std::to_string(NLi * Ntests + i + 1), outputFolder,
                                                     "g");
                        Info << "debug = " << NLi * noiseLevelVec.size() + i + 1 << endl;
                        Info << "************************************" << endl;
                        Info << endl << endl;
                    }
                    else
                    {
                        Info << "CG did not converged" << endl;
                        Info << "************************************" << endl;
                        i--;
                        Info << endl << endl;
                    }
                }
                else if(ParamBCnoiseLevelTest)
                {
                    example_paramBC.parameterizedBC("fullPivLU", 3);
                    volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
                    ITHACAstream::exportSolution(gParametrizedField,
                                                 std::to_string(NLi * Ntests + i + 1),
                                                 outputFolder,
                                                 "g");
                }
                example_paramBC.Tmeas = TmeasOrig;
                example_CG.Tmeas = TmeasOrig;
            }
        }

        example_paramBC.postProcess(outputFolder, "g");
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (CGandParamErrorTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing parameterized BC inverse solver with NOISY data" <<
             endl;
        word outputFolder = "./ITHACAoutput/CGandParamErrorTest/";
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<List<scalar>> heatFluxWeights;
        Eigen::VectorXd residualNorms;
        scalar innerField = 1.0;
        example_paramBC.set_gParametrized("rbf", 0.1);
        example_paramBC.parameterizedBCoffline();

        Info << "Introducing error in the measurements" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;

        
        Eigen::VectorXd varianceVect = example_paramBC.Tmeas * noiseLevel;
        varianceVect.cwiseProduct(varianceVect);
        example_CG.Jtol = example_CG.Tmeas.size() * varianceVect.maxCoeff() *
                       varianceVect.maxCoeff() / 2;
        std::cout << "Jtol = " << example_CG.Jtol << std::endl;
        example_CG.saveSolInLists = 1;

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example_paramBC.addNoise(noiseLevel);
            example_CG.addNoise(noiseLevel);

            //TSVD solution
            example_paramBC.parameterizedBC("TSVD", 3);
            volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(i * 3 + 1),
                                         outputFolder,
                                         "g");
            //LU solution
            example_paramBC.parameterizedBC("fullPivLU");
            gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(i * 3 + 2),
                                         outputFolder,
                                         "g");

            //CG solution
            example_CG.conjugateGradient();
            volScalarField heatFluxField = example_CG.list2Field(
                                               example_CG.gList[example_CG.gList.size() - 1]);
            ITHACAstream::exportSolution(heatFluxField,
                                         std::to_string(i * 3 + 3), outputFolder,
                                         "g");

            example_paramBC.Tmeas = TmeasOrig;
            example_CG.Tmeas = TmeasOrig;
        }

        example_paramBC.postProcess(outputFolder, "g", innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (parameterizedBCerrorTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing parameterized BC inverse solver with NOISY data" <<
             endl;
        word outputFolder = "./ITHACAoutput/parameterizedBCnoiseTest/";
        volScalarField gTrueField = example_paramBC.list2Field(example_paramBC.gTrue);
        ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");
        List<List<scalar>> heatFluxWeights;
        Eigen::VectorXd residualNorms;
        scalar innerField = 1.0;
        example_paramBC.set_gParametrized("rbf", 0.1);
        example_paramBC.parameterizedBCoffline();
        List<word> linSys_solvers;
        linSys_solvers.resize(1);
        //linSys_solvers[0] = "fullPivLU";
        //linSys_solvers[1] = "jacobiSvd";
        //linSys_solvers[2] = "householderQr";
        //linSys_solvers[3] = "ldlt";
        linSys_solvers[0] = "TSVD";
        Info << "Introducing error in the measurements" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        residualNorms.resize(Ntests * linSys_solvers.size());
        Eigen::VectorXd TmeasOrig = example_paramBC.Tmeas;
        auto density = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(1),Eigen::VectorXd::Ones(1));

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example_paramBC.addNoise(noiseLevel);

            List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
            List<scalar> solutionNorms;
            forAll(linSys_solvers, solverI)
            {
                Info << "Solver " << linSys_solvers[solverI] << endl;
                example_paramBC.parameterizedBC(linSys_solvers[solverI], 3);
                Info << endl;
                volScalarField gParametrizedField = example_paramBC.list2Field(example_paramBC.g);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
                                             outputFolder,
                                             "gParametrized");
                ITHACAstream::exportSolution(example_paramBC.T,
                                             std::to_string(i * linSys_solvers.size() + solverI + 1),
                                             outputFolder,
                                             "T");
                residualNorms(i * linSys_solvers.size() + solverI) = Foam::sqrt(
                            example_paramBC.residual.squaredNorm());
            }
            example_paramBC.Tmeas = TmeasOrig;
        }

        example_paramBC.postProcess(outputFolder, "gParametrized", innerField);
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

