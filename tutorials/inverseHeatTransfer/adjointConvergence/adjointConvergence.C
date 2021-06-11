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
    example_CG.of a heat transfer Reduction Problem
SourceFiles
    adjointConvergence.C
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
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    inverseLaplacianProblem_CG example_CG(argc, argv);
    //fvMesh& mesh(example_CG._mesh());
    //example_CG.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    //Setting parameters for the analytical benchmark
    double a = 5;
    double b = 10;
    double c = 15;
    double d = 20;
    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example_CG._mesh(),
                             example_CG._runTime());
    // Reading parameters from ITHACAdict
    example_CG.k =
        para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example_CG.k > 0, "thermalConductivity, k, not specified");
    example_CG.H =
        para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example_CG.H > 0, "Heat transfer coeff, H, not specified");
    example_CG.readThermocouples();
    example_CG.Tdiff = Eigen::VectorXd::Ones(example_CG.thermocouplesCellID.size())
                       * 55000;
    example_CG.solveAdjoint();
    volScalarField& lambda = example_CG._lambda();
    word folder = "./ITHACAoutput/lambda";
    ITHACAstream::exportSolution(lambda, "1", folder, "lambda");
    Info << "lambda L2norm = " << ITHACAutilities::L2Norm(lambda) << endl;
    return 0;
}

