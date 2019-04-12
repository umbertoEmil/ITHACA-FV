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
    14thermalBlock.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
// #include "laplacianProblem.H"
// #include "reducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>

// class tutorial14: public laplacianProblem
// {
// public:
// 	explicit tutorial14(int argc, char* argv[])
// 	:
// 	laplacianProblem(argc, argv),
// 	T(_T()),
// 	nu(_nu()),
// 	S(_S())
// 	{}
//         //! [tutorial02]
//         /// Temperature field
// 	volScalarField& T;
//         /// Diffusivity field
// 	volScalarField& nu;
//         /// Source term field
// 	volScalarField& S;
// };
// 

// void solve(volScalarField& T)
// {
// 	while (simple.loop())
//     {
//         Info<< "Time = " << runTime.timeName() << nl << endl;

//         while (simple.correctNonOrthogonal())
//         {
//             fvScalarMatrix TEqn
//             (
//                 fvm::ddt(T) - fvm::laplacian(DT, T)
//              ==
//                 fvOptions(T)
//             );

//             fvOptions.constrain(TEqn);
//             TEqn.solve();
//             fvOptions.correct(T);
//         }

//         #include "write.H"

//         runTime.printExecutionTime(Info);
//     }
// }


int main(int argc, char *argv[])
{
	// Optimization Variables
	int n_iter = 0;
	int n_iter_max = 10;

    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    volScalarField T0("T0",T);
    volScalarField f("T0",T);

    scalar v = 10.0;
    ITHACAutilities::assignIF(f,v);
    dimensionedScalar a("a", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1);





    // for (int i = 0; T.boundaryField().size(); i++)
    // {
    // 	if(T.boundaryField()[i].name() == "hotside")
    // 	{

    // 	}
    // }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;    

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::laplacian(DT, T) - f*a
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        // #include "write.H"
        // runTime.printExecutionTime(Info);
    }

    ITHACAstream::exportSolution(T, "2", "./results","Tout");


    ITHACAutilities::assignIF(T,v);
    ITHACAutilities::assignBC(T,2,v);

    surfaceScalarField heatFlux = fvc::snGrad(T);

    Eigen::MatrixXd Matrice;
    Eigen::VectorXd Vettore;
    List<Eigen::VectorXd> BCeigen;

    Vettore = Foam2Eigen::field2Eigen(T);
    BCeigen = Foam2Eigen::field2EigenBC(T);

    for (int i = 0; i<BCeigen.size(); i++)
    {
    std::cout << BCeigen[i] << std::endl;
	}

	// oggetto openfoam Info << OFOAM field << endl;
	// oggetto Eigen std::cout << eigen << std::endl;
	
	// export di un campo
	ITHACAstream::exportSolution(T, "1", "./results","Tout");

    // Info << T << endl;
    // Info << T0 << endl;    
    // ITHACAutilities::assignBC()


    // while(n_iter<n_iter_max)
    // {

    // }

    Info<< "End\n" << endl;

    return 0;
}
