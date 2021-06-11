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
    analyticalBenchmark_paramBCregularization.C
\*---------------------------------------------------------------------------*/

#include "MiniInverseProblems.h"
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
#include "../analyticalBenchmark.H"

using namespace SPLINTER;

class analyticalBenchmark_reg: public analyticalBenchmark
{
    public:
        explicit analyticalBenchmark_reg(int argc, char* argv[])
            :
            analyticalBenchmark(argc, argv)
        {
        }

        linearSystem linSys;
        Eigen::VectorXd Tmeas_err;
        Eigen::VectorXd sigma;



        void measError(double err_lev)
        {
            Eigen::VectorXd eta(Tmeas.size());
            double Tmeas_mean = Tmeas.sum() / Tmeas.size();
            sigma = err_lev * Tmeas;

            for (int i = 0; i < eta.size(); i++)
            {
                eta(i) = sigma(i) * stochastic::set_normal_random(0.0, 1.0);
            }

            Tmeas_err = Tmeas + eta;
            Info << "Variance max = " << sigma.maxCoeff() << endl;
            Info << "Variance min = " << sigma.minCoeff() << endl;
        }

        void getLinearSystem(Eigen::MatrixXd& A, Eigen::MatrixXd b)
        {
            linSys.A = A;
            linSys.b = b;
        }

        void LUpiv()
        {
            Info << "Usign LU with full pivoting" << endl;
            linSys.x = linSys.A.fullPivLu().solve(linSys.b);
            updateBC();
        }

        void TSVD(string regPar, double sigma, label truncatedSV = 0)
        {
            linSys.TSVD(regPar, sigma, truncatedSV);
            updateBC();
        }

        void tikhonov(string regPar, double sigma = 0.0, double alpha = 0.0)
        {
            linSys.tikhonov(regPar, sigma, alpha);
            updateBC();
        }

        void PCG()
        {
            Eigen::ConjugateGradient < Eigen::MatrixXd, Eigen::Lower | Eigen::Upper > cg;
            cg.setMaxIterations(2);
            cg.compute(linSys.A);
            linSys.x = cg.solve(linSys.b);
            std::cout << "#iterations:     " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;
            updateBC();
        }

        void updateBC()
        {
            gWeights.resize(linSys.x.rows());
            forAll(gWeights, weightI)
            {
                gWeights[weightI] = linSys.x(weightI);
            }
            update_gParametrized(gWeights);
            solveDirect();
        }
};

int main(int argc, char* argv[])
{
    solverPerformance::debug = 0; //No verbose output
    double time;
    analyticalBenchmark_reg example(argc, argv);
    //Setting parameters for the analytical benchmark
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    double W = 1; //Domain width [m]
    // Reading tests to perform
    ITHACAparameters para;
    label regularizationSetupTest =
        para.ITHACAdict->lookupOrDefault<int>("regularizationSetupTest", 0);
    label parameterizedBCerrorTest_TSVD =
        para.ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest_TSVD", 0);
    // Reading parameters from ITHACAdict
    example.thermocouplesNum =
        para.ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    M_Assert(example.thermocouplesNum > 0, "Number of thermocouples not specified");
    example.k = para.ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para.ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    label Ntests = para.ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double refGrad = para.ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction = para.ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);
    label TSVDregularization =
        para.ITHACAdict->lookupOrDefault<int>("TSVDregularization", 0);
    label TikhonovRegularization =
        para.ITHACAdict->lookupOrDefault<int>("TikhonovRegularization", 0);
    label PCGregularization =
        para.ITHACAdict->lookupOrDefault<int>("PCGregularization", 0);
    // Setting BC at the cold side
    example.coldSide_ind = example.mesh.boundaryMesh().findPatchID("coldSide");
    label coldSideSize = example.T.boundaryField()[example.coldSide_ind].size();
    example.Tf.resize(coldSideSize);
    example.refGrad.resize(coldSideSize);
    example.valueFraction.resize(coldSideSize);
    forAll(example.Tf, faceI)
    {
        scalar faceZ =
            example.mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].z();
        scalar faceX =
            example.mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].x();
        example.Tf[faceI] = example.k / example.H * b * faceX + a * faceX * faceX +
                            b * faceX * W - a * faceZ * faceZ + c;
        example.refGrad[faceI] = refGrad;
        example.valueFraction[faceI] = valueFraction;
    }
    // Setting BC at hotSide
    example.hotSide_ind = example.mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = example.T.boundaryField()[example.hotSide_ind].size();
    example.heatFlux_hotSide.resize(hotSideSize);
    forAll(example.heatFlux_hotSide, faceI)
    {
        scalar faceX =
            example.mesh.boundaryMesh()[example.hotSide_ind].faceCentres()[faceI].x();
        example.heatFlux_hotSide[faceI] = - b * faceX;
    }
    // Setting BC at gammaEx1
    example.gammaEx1_ind = example.mesh.boundaryMesh().findPatchID("gammaEx1");
    label gammaEx1Size = example.T.boundaryField()[example.gammaEx1_ind].size();
    example.heatFlux_gammaEx1.resize(gammaEx1Size);
    forAll(example.heatFlux_gammaEx1, faceI)
    {
        scalar faceZ =
            example.mesh.boundaryMesh()[example.gammaEx1_ind].faceCentres()[faceI].z();
        example.heatFlux_gammaEx1[faceI] = - 2 * a * faceZ;
    }
    // Setting BC at gammaEx2
    example.gammaEx2_ind = example.mesh.boundaryMesh().findPatchID("gammaEx2");
    label gammaEx2Size = example.T.boundaryField()[example.gammaEx2_ind].size();
    example.heatFlux_gammaEx2.resize(gammaEx2Size);
    forAll(example.heatFlux_gammaEx2, faceI)
    {
        scalar faceX =
            example.mesh.boundaryMesh()[example.gammaEx2_ind].faceCentres()[faceI].x();
        scalar faceY =
            example.mesh.boundaryMesh()[example.gammaEx2_ind].faceCentres()[faceI].y();
        example.heatFlux_gammaEx2[faceI] = (2 * a * faceX + b * faceY);
    }
    // Setting BC at gammaEx3
    example.gammaEx3_ind = example.mesh.boundaryMesh().findPatchID("gammaEx3");
    label gammaEx3Size = example.T.boundaryField()[example.gammaEx3_ind].size();
    example.heatFlux_gammaEx3.resize(gammaEx3Size);
    forAll(example.heatFlux_gammaEx3, faceI)
    {
        scalar faceZ =
            example.mesh.boundaryMesh()[example.gammaEx3_ind].faceCentres()[faceI].z();
        example.heatFlux_gammaEx3[faceI] =  (2 * a * faceZ);
    }
    // Setting BC at gammaEx4
    example.gammaEx4_ind = example.mesh.boundaryMesh().findPatchID("gammaEx4");
    label gammaEx4Size = example.T.boundaryField()[example.gammaEx4_ind].size();
    example.heatFlux_gammaEx4.resize(gammaEx4Size);
    forAll(example.heatFlux_gammaEx4, faceI)
    {
        scalar faceX =
            example.mesh.boundaryMesh()[example.gammaEx4_ind].faceCentres()[faceI].x();
        scalar faceY =
            example.mesh.boundaryMesh()[example.gammaEx4_ind].faceCentres()[faceI].y();
        example.heatFlux_gammaEx4[faceI] = - example.k * (2 * a * faceX + b * faceY);
    }
    example.solveTrue();
    // setting analytical solution
    volScalarField T_true(example.T);

    for (label i = 0; i < T_true.internalField().size(); i++)
    {
        auto cx = T_true.mesh().C()[i].component(vector::X);
        auto cy = T_true.mesh().C()[i].component(vector::Y);
        auto cz = T_true.mesh().C()[i].component(vector::Z);
        T_true.ref()[i] = a * cx * cx + b * cx * cy - a * cz * cz + c;
    }

    Info << "Exporting analytical solution" << endl;
    ITHACAstream::exportSolution(T_true, "1", "./ITHACAoutput/true/",
                                 "analyticalSol");
    // Setting up the thermocouples
    example.readThermocouples();
    forAll(example.thermocouplesCellID, cellI)
    {
        example.Tmeas(cellI) =
            example.T.internalField()[example.thermocouplesCellID[cellI]];
    }

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

    if (regularizationSetupTest)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing the implementation of regularizing techniques for parametrized BC inverse solver"
             <<
             endl;
        word folder = "./ITHACAoutput/regularizationSetupTest/";
        double residualNorms;
        scalar innerField = 0.0;
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        Eigen::VectorXd variance = example.Tmeas * 0.02;
        Info << "Noise variance mean = " << variance.mean() << endl;
        Eigen::VectorXd TmeasOrig = example.Tmeas;
        example.Tmeas = TmeasOrig;
        Eigen::VectorXd measurementsError(example.Tmeas.size());

        for (int i = 0; i < example.Tmeas.size(); i++)
        {
            measurementsError(i) = example.Tmeas.mean() * 0.02 *
                                   stochastic::set_normal_random(0.0, 1.0);
        }

        example.Tmeas += measurementsError;
        Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
        Eigen::MatrixXd b = example.Theta.transpose() * (example.Tmeas +
                            example.addSol);
        example.getLinearSystem(A, b);
        Info << "example.gWeights.size() = " << example.gWeights.size() << endl;
        std::string regPar = "UPRE";
        // LU with pivoting
        example.LUpiv();
        volScalarField gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(1),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(1),
                                     folder,
                                     "T");
        // TSVD
        example.TSVD(regPar, measurementsError.mean(), 3);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(2),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(2),
                                     folder,
                                     "T");
        // Tikhonov
        example.tikhonov(regPar, measurementsError.maxCoeff());
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(3),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(3),
                                     folder,
                                     "T");
        // Tikhonov
        example.tikhonov("GCV", measurementsError.maxCoeff());
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(4),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(4),
                                     folder,
                                     "T");
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(example.linSys.A,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd singularValues = svd.singularValues();
        ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen",
                                   folder);
        std::cout << "Residual 2-norm = " << residualNorms << std::endl;
        example.postProcess(folder, "gParametrized");
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (parameterizedBCerrorTest_TSVD)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing TSVD parameterized BC inverse solver with NOISY data" <<
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
        Info << "Introducing error in the measurements" << endl;
        Info << "Performing " << Ntests << " tests." << endl;
        residualNorms.resize(Ntests);
        Eigen::VectorXd TmeasOrig = example.Tmeas;

        for (label i = 0; i < Ntests; i++)
        {
            Info << "Test " << i << endl;
            example.Tmeas = TmeasOrig;
            Eigen::VectorXd measurementsError(example.Tmeas.size());

            for (int i = 0; i < example.Tmeas.size(); i++)
            {
                measurementsError(i) = example.Tmeas.mean() * 0.02 *
                                       stochastic::set_normal_random(0.0, 1.0);
            }

            example.Tmeas += measurementsError;
            Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
            Eigen::MatrixXd b = example.Theta.transpose() * (example.Tmeas +
                                example.addSol);
            example.getLinearSystem(A, b);
            List<List<scalar>> heatFluxWeights_err = heatFluxWeights;
            List<scalar> solutionNorms;

            if (TSVDregularization)
            {
                Info << "Using TSVD" << endl;
                example.parameterizedBC(outputFolder, "TSVD", 2);
            }
            else if (TikhonovRegularization)
            {
                Info << "Using Tikhonov" << endl;
                example.tikhonov("UPRE", measurementsError.maxCoeff(), 0.00001);
            }
            else if (PCGregularization)
            {
                Info << "Using PCG" << endl;
                example.PCG();
            }

            volScalarField gParametrizedField = example.list2Field(example.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "gParametrized");
            ITHACAstream::exportSolution(example.T,
                                         std::to_string(i + 1),
                                         outputFolder,
                                         "T");
            residualNorms(i) = Foam::sqrt(
                                   example.residual.squaredNorm());
            Info << "Measurements error L2 norm= " << measurementsError.norm() <<
                 endl;
        }

        example.parameterizedBC_postProcess(outputFolder, innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    return 0;
}

