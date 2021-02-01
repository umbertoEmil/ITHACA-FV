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
    analyticalBenchmark_unsteady.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "pimpleControl.H"
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

#include "unsteadyTest.H"
using namespace SPLINTER;


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    unsteadyTest example(argc, argv);
    // Reading tests to perform
    ITHACAparameters para;
    example.k = para.ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.rho = para.ITHACAdict->lookupOrDefault<double>("density", 0);
    M_Assert(example.rho > 0, "Density, rho, not specified");
    example.Cp = para.ITHACAdict->lookupOrDefault<double>("heatCapacity", 0);
    M_Assert(example.Cp > 0, "heatCapacity, Cp, not specified");
    example.H = para.ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    example.alpha = para.ITHACAdict->lookupOrDefault<double>("diffusivity", 0);
    M_Assert(example.alpha > 0, "Diffusivity, alpha, not specified");
    double refGrad = para.ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction = para.ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);
    fvMesh& mesh = example._mesh();
    volScalarField& T = example._T();
    

    //example.readThermocouples();
    
    // Setting BC at the cold side
    example.coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
    label coldSideSize = T.boundaryField()[example.coldSide_ind].size();
    example.Tf.resize(coldSideSize);
    example.refGrad.resize(coldSideSize);
    example.valueFraction.resize(coldSideSize);
    forAll(example.Tf, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].x();
        scalar faceY =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].y();
        scalar faceZ =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].z();

        example.Tf[faceI] = 20; 
        example.refGrad[faceI] = refGrad;
        example.valueFraction[faceI] = valueFraction;
    }
    // Setting BC at hotSide
    example.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T.boundaryField()[example.hotSide_ind].size();
    example.heatFlux_hotSide.resize(hotSideSize);
    forAll(example.heatFlux_hotSide, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].x();
        scalar faceY =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].y();
        scalar faceZ =
            mesh.boundaryMesh()[example.coldSide_ind].faceCentres()[faceI].z();
        example.heatFlux_hotSide[faceI] = 100;
    }

    // Performing unsteady  solution
    auto t1 = std::chrono::high_resolution_clock::now();
    example.solveUnsteady();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6;
    std::cout << "Unsteady solution took  = " << duration << " seconds" << std::endl;
    
    return 0;
}

