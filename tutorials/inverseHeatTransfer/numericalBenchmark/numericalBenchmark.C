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




int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    double time;
    numericalBenchmark_CG example_CG(argc, argv);
    numericalBenchmark_paramBC example_paramBC(argc, argv);
    // Reading parameters from ITHACAdict
    ITHACAparameters* para = ITHACAparameters::getInstance(example_paramBC._mesh(),
                             example_paramBC._runTime());
    example_CG.g_0 = para->ITHACAdict->lookupOrDefault<double>("g_0", 0.0);
    example_CG.g_X = para->ITHACAdict->lookupOrDefault<double>("g_X", 0.0);
    example_CG.g_Z = para->ITHACAdict->lookupOrDefault<double>("g_Z", 0.0);
    example_CG.Tf_0 = para->ITHACAdict->lookupOrDefault<double>("Tf_0", 0.0);
    example_CG.Tf_delta = para->ITHACAdict->lookupOrDefault<double>("Tf_delta",
                          0.0);
    example_paramBC.g_0 = example_CG.g_0;
    example_paramBC.g_X = example_CG.g_X;
    example_paramBC.g_Z = example_CG.g_Z;
    example_paramBC.Tf_0 = example_CG.Tf_0;
    example_paramBC.Tf_delta = example_CG.Tf_delta;
    word solver = para->ITHACAdict->lookupOrDefault<word>("linSysSolver",
                  "fullPivLU");
    label TSVDtruncation =
        para->ITHACAdict->lookupOrDefault<label>("TSVDtruncation", 3);
    label CGtest = para->ITHACAdict->lookupOrDefault<int>("CGtest", 0);
    label CGnoiseTest = para->ITHACAdict->lookupOrDefault<int>("CGnoiseTest", 0);
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    label parameterizedBCtest_RBFwidth =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest_RBFwidth", 0);
    label parameterizedBC_noiseLevelTest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBC_noiseLevelTest", 0);
    label parameterizedBCerrorTest_TSVD =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCerrorTest_TSVD", 0);
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
    double RBFshapePar = para->ITHACAdict->lookupOrDefault<double>("RBFshapePar",
                         0);
    int TSVDtrunc = para->ITHACAdict->lookupOrDefault<int>("TSVDregularization", 0);
    int rbfWidthTest_size =
        para->ITHACAdict->lookupOrDefault<int>("rbfWidthTest_size", 0);
    double noiseLevel = para->ITHACAdict->lookupOrDefault<double>("noiseLevel",
                        0);
    example_CG.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity",
                   0);
    M_Assert(example_CG.k > 0, "thermalConductivity, k, not specified");
    example_CG.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    example_paramBC.k = example_CG.k;
    example_paramBC.H = example_CG.H;
    M_Assert(example_CG.H > 0, "Heat transfer coeff, H, not specified");
    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double refGrad = para->ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction =
        para->ITHACAdict->lookupOrDefault<double>("valueFraction",
                0.0);
    //*******************************************//
    fvMesh& mesh = example_paramBC._mesh();
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
        example_CG.g[faceI] = (example_CG.g_X * (faceX - 1) * (faceX - 1) +
                               example_CG.g_Z * faceZ + example_CG.g_0) ;
    }
    example_CG.set_Tf();
    example_CG.gTrue = example_CG.g;
    example_CG.solveTrue();
    example_paramBC.g = example_CG.g;
    example_paramBC.gTrue = example_CG.gTrue;
    example_paramBC.set_Tf();
    example_paramBC.solveTrue();
    volScalarField& T(example_CG._T());
    // Setting up the thermocouples
    example_CG.readThermocouples();
    example_CG.Tmeas = example_CG.fieldValueAtThermocouples(T);
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
#include"CGtest.H"
    }

    if (parameterizedBCtest)
    {
#include"parameterizedBCtest.H"
    }

    if (parameterizedBCtest_RBFwidth)
    {
#include"parameterizedBCtest_RBFwidth.H"
    }

    if (parameterizedBC_noiseLevelTest)
    {
#include"paramBC_noiseLevelTest.H"
    }

    if (parameterizedBCerrorTest_TSVD)
    {
#include"parameterizedBCerrorTest_TSVD.H"
    }

    return 0;
}


