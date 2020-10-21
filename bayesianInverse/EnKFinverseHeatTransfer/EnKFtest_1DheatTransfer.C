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
    EnKFtest_1DheatTransfer.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"

#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "muq2ithaca.H"

#include "EnKFtest_1DheatTransfer.H"
using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    int Nseeds = 50;
    EnKFtest_1DheatTransfer example(argc, argv, Nseeds);
    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                         example._runTime());

    example.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.rho = para->ITHACAdict->lookupOrDefault<double>("density", 0);
    M_Assert(example.rho > 0, "Density, rho, not specified");
    example.Cp = para->ITHACAdict->lookupOrDefault<double>("heatCapacity", 0);
    M_Assert(example.Cp > 0, "heatCapacity, Cp, not specified");

    fvMesh& mesh = example._mesh();
    volScalarField& T = example._T();
    int stateSize = T.size();
    scalar initialField = 10; //TODO make it coherent with dictionary value

        
    example.solveDirect();
    int Ntimes = example.Ntimes;

    // Setting up the densities
    // Gaussian densities are assumed
    example.priorSetup(initialField, 0.7);
    example.modelErrorSetup(0.0, 0.7);
    example.measNoiseSetup(0.0, 0.05);

    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);

    example.priorSampling();


    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    Eigen::MatrixXd minConfidence = posteriorMean;
    Eigen::MatrixXd maxConfidence = minConfidence;

    posteriorMean.col(0) = posteriorSamples.rowwise().mean();

    example.reconstruction();
    
    return 0;
}

