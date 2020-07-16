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
    testFieldReconstruction.C
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
#include "testFieldReconstruction.H"

// The objective of this code is to determine the heat flux "g" on the boundary "hotSide"

int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    double time;
    testFieldReconstruction example(argc, argv);
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
    label parameterizedBCtest =
        para->ITHACAdict->lookupOrDefault<int>("parameterizedBCtest", 0);
    example.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");

    word outputFolder = "./ITHACAoutput/parameterizedBCtest/";

    // Setting up the true problem with known heat flux g
    fvMesh& mesh = example._mesh();
    volScalarField& T(example._T());
    example.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T.boundaryField()[example.hotSide_ind].size();
    example.g.resize(hotSideSize);
    example.gTrue.resize(hotSideSize);
    // Here I define the heat flux that I want to estimate
    forAll(example.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example.hotSide_ind].faceCentres()[faceI].x();
        example.g[faceI] = example.k * (b * faceX + c) ;
    }
    example.gTrue = example.g; 
    example.solveTrue(outputFolder);
    

    // Setting up the thermocouples and filling the vector containing the 
    // measurements Tmeas with the values of the temperature field computed
    // by solveTrue 
    example.readThermocouples();
    example.Tmeas = example.fieldValueAtThermocouples(T);

    // Solving inverse problem
    Info << endl;
    Info << "*********************************************************" << endl;
    Info << "Performing test for the parameterized BC inverse solver" << endl;
    Info << endl;
    volScalarField gTrueField = example.list2Field(example.gTrue);
    ITHACAstream::exportSolution(gTrueField,
                                 "1", outputFolder,
                                 "gTrue");

    
    // Now I compute the basis  
    scalar rbfShapeParameter = 0.7;
    example.set_gParametrized("rbf", rbfShapeParameter); // define the basis of the heat flux
    label forceComputingBasis = 1;
    example.parameterizedBCoffline(forceComputingBasis); 
    PtrList<volScalarField> Tad_base1 = example.Tad_base;
    PtrList<volScalarField> Tbasis1 = example.Tbasis;
    Eigen::VectorXd weightsTry(16);
    
    //These are weights that I made up
    //weightsTry << -10000, +15000, 0, -30000, 8000, -15000, 0, 45000, 0, 15000, 0, 0, -8000, 15000, 0, -40500;

    //These are weights that the inverse solver gives as a result
    weightsTry << -14460.3, 21096.3, 0, -6924.19, 21020.8, -31081.5, 0, 10093.6, 0, 814.007, 0, 0, -7179.26, 10049.4, 0, -3427.73;

    // Reconstruct the field using weightsTry
    example.parameterizedBCpostProcess(weightsTry);
    volScalarField& Trecon1(example._T());
    ITHACAstream::exportSolution(Trecon1,
                     std::to_string(1),
                     outputFolder,
                     "Trecon1");

    // Now I read the basis from file and reconstruct the field
    Info << " \n\n AGAIN \n\n";
    example.parameterizedBCoffline(); 
    PtrList<volScalarField> Tad_base2 = example.Tad_base;
    PtrList<volScalarField> Tbasis2 = example.Tbasis;
    example.parameterizedBCpostProcess(weightsTry);
    volScalarField& Trecon2(example._T());
    ITHACAstream::exportSolution(Trecon2,
                     std::to_string(1),
                     outputFolder,
                     "Trecon2");

    // Compute the difference between computed and read basis and write it
    forAll(Tbasis1, baseI)
    {
        volScalarField test = Tbasis1[baseI] - Tbasis2[baseI];
	ITHACAstream::exportSolution(test,
                 std::to_string(baseI),
                 outputFolder,
                 "differenceBetweenBasis");

    }


    //------------ description of offline phase -------------//
    // Let phi_i be the basis of the heat flux g
    // In the offline phase I whant to compute T[phi_i]
    // for all i. Moreover, I have to compute what I call 
    // the additional problem, Tadd, wich do not depend on phi_i.
    // Then, I have to take the values of T[phi_i] and Tadd at 
    // the termocouples and store them in vectors. With these 
    // values, I assemble the matrix Theta which is such that
    //              
    //     Theta[i,j] = T[phi_j](thermocouple_i) 
    //                  + Tadd(thermocouple_i)  
    //
    // Now I have everything I need for the online phase
    // 

    // Online
    //------------ description of online phase -------------//
    // In the online phase, I want to determine the weights
    // of the basis function that minimize the cost function
    //
    //     J = sum_i (Tcomp_i - Tmeas_I)^2
    //
    // thanks to the affinity of the tmeperature field with 
    // respect to the basis of the heat flux, I can find 
    // these weights by solving
    //
    //     Theta^T Theta gWeights = Theta^T (Tmeas + Tadd)
    //
    // In parameterizedBC, I solve this linear system and then
    // I reconstruct the respective temperature field by
    //
    //     T[gWeights] = sum_i gWeights_i (T[phi_i] + Tadd)
    //                   - Tadd
    //
    // this is done by the function reconstructT 
    //
    //label dummy = 0;
    //example.parameterizedBC("fullPivLU", dummy); 

    //label dummy = 0;
    //Eigen::VectorXd weights = example.parameterizedBC("fullPivLU", dummy);
    //
    //std::cout << "\n weights = \n " << weights.transpose() << std::endl;
    //std::cout << "\n weights diff = \n " << (weightsTry - weights).transpose() << std::endl;


    //// Results output
    //volScalarField gParametrizedField = example.list2Field(example.g);
    //ITHACAstream::exportSolution(gParametrizedField,
    //                             std::to_string(1),
    //                             outputFolder,
    //                             "gParametrized");
    //volScalarField& Tout(example._T());
    //ITHACAstream::exportSolution(Tout,
    //                             std::to_string(1),
    //                             outputFolder,
    //                             "T");

    Info << "*********************************************************" << endl;
    Info << endl;

    return 0;
}

