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
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
// #include "reducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"

using namespace SPLINTER;


class Tutorial15: public inverseLaplacianProblem
{
    public:
        explicit Tutorial15(int argc, char* argv[])
            :
            inverseLaplacianProblem(argc, argv),
            T(_T()),
            lambda(_lambda()),
            deltaT(_deltaT()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }
        volScalarField& T;
        volScalarField& lambda;
        volScalarField& deltaT;
        fvMesh& mesh;
        Time& runTime;

        List<scalar> gTrue;

        void solveTrue()
        {
            g.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
            gTrue.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
            forAll (T.boundaryField()[hotSide_ind], faceI)
            {
                gTrue[faceI] = - 100000 *
                               mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z() *
                               mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
            }
            valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
            valueFractionAdj.resize(mesh.boundaryMesh()["coldSide"].size());
            Eigen::VectorXd faceCellDist =
                ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
            forAll (valueFraction, faceI)
            {
                valueFraction[faceI] =  1 / (1 + (k / H / faceCellDist(faceI)));
                valueFractionAdj[faceI] =  1 / (1 + (1 / k / H / faceCellDist(faceI)));
            }
            // True problem BCrefGrad
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(T, patchI, - gTrue / k);
                }
                else
                {
                    ITHACAutilities::assignBC(T, patchI, homogeneousBC);
                }
            }
            fvMesh& mesh = _mesh();
            volScalarField& T = _T();
            Foam::Time& runTime = _runTime();
            ITHACAutilities::assignIF(T, homogeneousBC);
            simpleControl simple(mesh);
#if OFVER == 6

            while (simple.loop(runTime))
#else
            while (simple.loop())
#endif
            {
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(DT, T)
                    );
                    TEqn.solve();
                }
            }

            //Reinitializing runTime
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 1);
            //solve("direct");
            readThermocouples();
            forAll(thermocouplesCellID, cellI)
            {
                Tmeas(cellI) = T.internalField()[thermocouplesCellID[cellI]];
            }

            if (interpolation)
            {
                Info << "Interpolating thermocouples measurements in the " <<
                     "plane defined by the thermocouples" << endl;
                thermocouplesInterpolation();
            }
            else
            {
                Info << "NOT interpolating thermocouples measurements" << endl;
            }
        }

};






int main(int argc, char* argv[])
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double time;
    Tutorial15 example(argc, argv);
    ITHACAparameters para;
    example.cgIterMax = para.ITHACAdict->lookupOrDefault<int>("cgIterMax", 100);
    example.thermocouplesNum =
        para.ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    M_Assert(example.thermocouplesNum > 0, "Number of thermocouples not specified");
    example.interpolation = para.ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example.Jtol =  para.ITHACAdict->lookupOrDefault<double>("Jtolerance",
                    0.000001);
    example.JtolRel =
        para.ITHACAdict->lookupOrDefault<double>("JrelativeTolerance",
                0.001);
    example.k = para.ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para.ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    double Tf = para.ITHACAdict->lookupOrDefault<double>("Tf", 300.0);
    double refGrad = para.ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double homogeneousBClist = 0.0;
    double valueFraction = para.ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);
    example.Tf.resize(example.T.boundaryField()[example.hotSide_ind].size());
    example.refGrad.resize(example.T.boundaryField()[example.hotSide_ind].size());
    example.valueFraction.resize(
        example.T.boundaryField()[example.hotSide_ind].size());
    forAll(example.Tf, faceI)
    {
        example.Tf[faceI] = Tf;
        example.refGrad[faceI] = refGrad;
        example.valueFraction[faceI] = valueFraction;
    }
    //example.solveTrue();
    //ITHACAstream::exportSolution(example.T, "2", "./ITHACAoutput/Offline/", "Ttrue");
    example.readThermocouples();
    Eigen::MatrixXd D;
    cnpy::load(D, "./thermocouplesSamples.npy");
    example.Tmeas.resize(D.rows());
    example.Tmeas = D.col(0);
    example.thermocouplesInterpolation();

    if (example.conjugateGradient())
    {
        Info << "Exiting conjugate gradient" << endl;
    }
    else
    {
        Info << "Max iterations condition met\n" << endl;
    }

    example.writeFields(1, "./results");
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>
                     (t2 - t1);
    time = time_span.count();
    std::cout << "CPU time = " << time << std::endl;
    reduce(time, maxOp<double>());
    return 0;
}
