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
#include <float.h>
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
    label CGnoiseLevelTest =
        para->ITHACAdict->lookupOrDefault<int>("CGnoiseLevelTest", 0);
    label ParamBCnoiseLevelTest =
        para->ITHACAdict->lookupOrDefault<int>("ParamBCnoiseLevelTest", 0);
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
    label thermocouplesLocationTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_CG", 0);
    label thermocouplesLocationTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesLocationTest_paramBC", 0);
    label thermocouplesNumberTest_CG =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_CG", 0);
    label thermocouplesNumberTest_paramBC =
        para->ITHACAdict->lookupOrDefault<int>("thermocouplesNumberTest_paramBC", 0);
    // Reading parameters from ITHACAdict
    example_CG.cgIterMax = para->ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example_CG.interpolation =
        para->ITHACAdict->lookupOrDefault<int>("interpolation",
                1);
    example_CG.Jtol =  para->ITHACAdict->lookupOrDefault<double>("Jtolerance",
                       0.000001);
    example_CG.JtolRel =
        para->ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);
    double rbfShapePar = para->ITHACAdict->lookupOrDefault<double>("rbfShapePar",
                         0);
    M_Assert(rbfShapePar > 0, "rbfShapePar not specified");
    example_paramBC.k =
        para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_paramBC.k > 0, "thermalConductivity, k, not specified");
    example_paramBC.H =
        para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_paramBC.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.k = example_paramBC.k;
    example_CG.H = example_paramBC.H;
    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double noiseLevel = para->ITHACAdict->lookupOrDefault<double>("noiseLevel",
                        0);
    int rbfWidthTest_size =
        para->ITHACAdict->lookupOrDefault<int>("rbfWidthTest_size", 0);
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
    Info << "L2 norm of the error = " << ITHACAutilities::errorL2Abs(T_true,
            T) << endl;
    Info << "Linf norm of the error = " << ITHACAutilities::LinfNorm(error) << endl;
    Info << "L2 norm of T = " << ITHACAutilities::L2Norm(T) << endl;
    Info << "L2 norm of Tanal = " << Ttrue_L2norm << endl;
    Info << "L2 norm of diff = " << ITHACAutilities::L2Norm(
             T) - Ttrue_L2norm  <<  endl;
    Info << "L2 norm of rel diff = " << (ITHACAutilities::L2Norm(
            T) - Ttrue_L2norm) / Ttrue_L2norm  << endl << endl;
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
#include "CGtest.H"
    }

    if (parameterizedBCtest)
    {
#include"parameterizedBCtest.H"
    }

    if (parameterizedBCtest_RBFwidth)
    {
#include"parameterizedBCtest_RBFwidth.H"
    }

    if (CGnoiseTest)
    {
#include"CGnoiseTest.H"
    }

    if (CGnoiseLevelTest || ParamBCnoiseLevelTest)
    {
#include"noiseLevelTest.H"
    }

    if (CGandParamErrorTest)
    {
#include"CGandParamErrorTest.H"
    }

    if (parameterizedBCerrorTest)
    {
#include"parameterizedBCerrorTest.H"
    }

    if (parameterizedBCerrorTest_TSVD)
    {
#include"parameterizedBCerrorTest_TSVD.H"
    }

    if (thermocouplesLocationTest_CG)
    {
#include"thermocouplesLocation_CG.H"
    }

    if (thermocouplesLocationTest_paramBC)
    {
#include"thermocouplesLocation_paramBC.H"
    }

    if (thermocouplesNumberTest_CG)
    {
#include"thermocouplesNumberTest_CG.H"
    }

    if (thermocouplesNumberTest_paramBC)
    {
#include"thermocouplesNumberTest_paramBC.H"
    }

    return 0;
}

