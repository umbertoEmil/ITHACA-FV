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
    18bayesian.C
\*---------------------------------------------------------------------------*/


//#include "ITHACAbayesian.H"
#include "MiniInverseProblems.h"
#include <iostream>
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"

using namespace SPLINTER;

class Tutorial18: public inverseLaplacianProblem
{
    public:
        explicit Tutorial18(int argc, char* argv[])
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

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> lambdaFullField;
        PtrList<volScalarField> deltaTfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gDiffField;
        PtrList<volScalarField> gRelErrField;

        linearSystem linSys;
        Eigen::VectorXd Tmeas_err;
        double sigma;

        scalar L2norm_g(volScalarField& field1, volScalarField& field2)
        {
            scalar L2 = 0;
            //Access the mesh information for the boundary
            const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
            //List of cells close to a boundary
            const labelUList& faceCells = cPatch.faceCells();
            forAll(cPatch, faceI)
            {
                scalar faceX = mesh.Cf().boundaryField()[hotSide_ind][faceI][0];
                scalar faceZ = mesh.Cf().boundaryField()[hotSide_ind][faceI][2];

                if (isInPlane(faceX, faceZ))
                {
                    //id of the owner cell having the face
                    label faceOwner = faceCells[faceI] ;
                    scalar faceArea = mesh.magSf().boundaryField()[hotSide_ind][faceI];
                    scalar relErr = (field1[faceOwner] - field2[faceOwner]) / field2[faceOwner];
                    L2 += faceArea * relErr * relErr;
                }
            }
            return Foam::sqrt(L2);
        }

        scalar LinfinityErrorNorm(volScalarField& field1, volScalarField& field2)
        {
            scalar Linfty = 0;
            //Access the mesh information for the boundary
            const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
            //List of cells close to a boundary
            const labelUList& faceCells = cPatch.faceCells();
            label firstElement = 1;
            forAll(cPatch, faceI)
            {
                scalar faceX = mesh.Cf().boundaryField()[hotSide_ind][faceI][0];
                scalar faceZ = mesh.Cf().boundaryField()[hotSide_ind][faceI][2];

                if (isInPlane(faceX, faceZ))
                {
                    //id of the owner cell having the face
                    label faceOwner = faceCells[faceI] ;
                    scalar fieldsDiff = field1[faceOwner] - field2[faceOwner];
                    scalar value = fieldsDiff * fieldsDiff / (field2[faceOwner] *
                                   field2[faceOwner]);
                    value = Foam::sqrt(value);

                    if (firstElement)
                    {
                        Linfty = value;
                        firstElement = 0;
                    }
                    else if (Linfty < value)
                    {
                        Linfty = value;
                    }
                }
            }
            return Linfty;
        }

        void computeRelativeErrorFields()
        {
            //Access the mesh information for the boundary
            const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
            //List of cells close to a boundary
            const labelUList& faceCells = cPatch.faceCells();
            gRelErrField = gField;
            forAll(gField, testI)
            {
                forAll(cPatch, faceI)
                {
                    label faceOwner = faceCells[faceI] ;
                    gRelErrField[testI][faceOwner] = (gField[testI][faceOwner] -
                                                      gFullField[0][faceOwner]) /  (gFullField[0][faceOwner]);
                }
            }
        }

        PtrList<volScalarField> differenceField (PtrList<volScalarField>& field1,
                PtrList<volScalarField>& field2)
        {
            PtrList<volScalarField> fieldDiff(field1);

            if (field2.size() == fieldDiff.size())
            {
                for (label i = 0; i < fieldDiff.size(); i++)
                {
                    fieldDiff[i] = (field1[i] - field2[i]);
                }
            }
            else if (field2.size() == 1)
            {
                for (label i = 0; i < fieldDiff.size(); i++)
                {
                    fieldDiff[i] = (field1[i] - field2[0]);
                }
            }
            else
            {
                Info << "Don't know how to behave in differenceField" << endl;
                Info << "the two fields have different size" << endl;
                Info << "Size field1 = " << field1.size()
                     << "Size field2 = " << field2.size() << endl;
            }

            return fieldDiff;
        }

        PtrList<volScalarField> relativeDifferenceField (PtrList<volScalarField>&
                field1, PtrList<volScalarField>& field2)
        {
            PtrList<volScalarField> fieldDiff;

            if (field2.size() == fieldDiff.size())
            {
                for (label i = 0; i < field1.size(); i++)
                {
                    fieldDiff.append((field1[i] - field2[i]) / field2[i]);
                }
            }
            else if (field2.size() == 1)
            {
                for (label i = 0; i < field1.size(); i++)
                {
                    fieldDiff.append((field1[i] - field2[0]) / field2[0]);
                }
            }
            else
            {
                Info << "WARNING:" << endl;
                Info << "Don't know how to behave in relativeDifferenceField" << endl;
                Info << "the two fields have different size" << endl;
                Info << "Size field1 = " << field1.size()
                     << ", Size field2 = " << field2.size() << endl;
            }

            return fieldDiff;
        }

        void computeTdiff()
        {
            TdiffField = TredField;
            forAll(TdiffField, testI)
            {
                TdiffField[testI] = TredField[testI] - TfullField[0];
            }
        }

        int isInPlane(double cx, double cz)
        {
            return (cx >= interpolationPlane.minX -
                    interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz >= interpolationPlane.minZ - interpolationPlane.thermocoupleCellDim[2] / 4 &&
                    cx <= interpolationPlane.maxX + interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz <= interpolationPlane.maxZ + interpolationPlane.thermocoupleCellDim[2] / 4
                   );
        }

        // To be added to ITHACAutilities
        scalar list2norm(List<scalar> list)
        {
            scalar norm = 0;
            forAll(list, I)
            {
                norm += list[I] * list[I];
            }
            return Foam::sqrt(norm);
        }

        void parameterizedBC_postProcess(word folder, scalar innerField = 0.0)
        {
            Info << "Computing errors" << endl;
            gFullField.resize(0);
            TfullField.resize(0);
            ITHACAstream::read_fields(mesh, gFullField, "gParametrized",
                                      folder);
            ITHACAstream::read_fields(mesh, TfullField, "T",
                                      folder);
            solveTrue();
            PtrList<volScalarField> gTrueField;
            gTrueField.resize(0);
            gTrueField.append(list2Field(gTrue, innerField));
            gDiffField = relativeDifferenceField(gFullField,
                                                 gTrueField);
            forAll(gDiffField, solutionI)
            {
                forAll(gDiffField[solutionI], faceI)
                {
                    gDiffField[solutionI][faceI] =
                        Foam::sqrt(gDiffField[solutionI][faceI] *
                                   gDiffField[solutionI][faceI]);
                }
            }
            Eigen::MatrixXd errorG_L2norm;
            errorG_L2norm.resize(gFullField.size(), 1);
            Eigen::MatrixXd errorG_LinfNorm = errorG_L2norm;
            gTrueField.resize(0);

            for (int i = 0; i < errorG_L2norm.rows() ; i++)
            {
                gTrueField.append(list2Field(gTrue, innerField));
                errorG_L2norm(i, 0) = L2norm_g(gFullField[i], gTrueField[0]);
                errorG_LinfNorm(i, 0) = LinfinityErrorNorm(gFullField[i],
                                        gTrueField[0]);
            }

            ITHACAstream::exportMatrix(errorG_L2norm, "relError_L2norm", "eigen",
                                       folder);
            ITHACAstream::exportMatrix(errorG_LinfNorm, "relError_LinfNorm", "eigen",
                                       folder);
            ITHACAstream::exportFields(gDiffField, folder, "gDiff");
            ITHACAstream::exportFields(gTrueField, folder, "gTrue");
        }

        void measError(double err_lev)
        {
            Eigen::VectorXd eta(Tmeas.size());
            double Tmeas_mean = Tmeas.sum() / Tmeas.size();
            sigma = err_lev / 100 * Tmeas_mean;

            for (int i = 0; i < eta.size(); i++)
            {
                eta(i) = sigma * stochastic::set_normal_random(0.0, 1.0);
            }

            Tmeas_err = Tmeas + eta;
        }

        void getLinearSystem(Eigen::MatrixXd& A, Eigen::MatrixXd b)
        {
            linSys.A = A;
            linSys.b = b;
        }

        void LUpiv()
        {
            linSys.x = linSys.A.fullPivLu().solve(linSys.b);
            updateBC();
        }

        void TSVD(string regPar, double sigma, label truncatedSV = 0)
        {
            linSys.TSVD(regPar, sigma, truncatedSV);
            updateBC();
        }

        void tikhonov(string regPar, double sigma = 0.0)
        {
            linSys.tikhonov(regPar, sigma);
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
    Tutorial18 example(argc, argv);
    // Reading parameters from ITHACAdict
    ITHACAparameters para;
    // Tests to do
    label fullOrderTest = para.ITHACAdict->lookupOrDefault<int>("fullOrderTest", 0);
    example.thermocouplesNum =
        para.ITHACAdict->lookupOrDefault<int>("thermocouplesNumber", 0);
    M_Assert(example.thermocouplesNum > 0, "Number of thermocouples not specified");
    // ADD read termocouples position from file
    example.interpolation = para.ITHACAdict->lookupOrDefault<int>("interpolation",
                            1);
    example.k = para.ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para.ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    double Tf = para.ITHACAdict->lookupOrDefault<double>("Tf", 300.0);
    double refGrad = para.ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction = para.ITHACAdict->lookupOrDefault<double>("valueFraction",
                           0.0);
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
        example.Tf[faceI] = Tf + 100 * Foam::sqrt(1.2 - faceZ);
        //example.Tf[faceI] = Tf;
        example.refGrad[faceI] = refGrad;
        example.valueFraction[faceI] = valueFraction;
    }

    if (1)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing the implementation of regularizing techniques for parametrized BC inverse solver"
             <<
             endl;
        word folder = "./ITHACAoutput/TSVD/";
        double residualNorms;
        scalar innerField = 1.0;
        double err_lev = 2; //% measurements error
        example.solveTrue();
        ITHACAstream::exportSolution(example.T, std::to_string(1),
                                     folder,
                                     "Ttrue");
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        example.measError(err_lev);
        Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
        Eigen::MatrixXd b = example.Theta.transpose() * (example.Tmeas_err +
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
        example.TSVD(regPar, example.sigma);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(2),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(2),
                                     folder,
                                     "T");
        // Tikhonov
        example.tikhonov(regPar, example.sigma);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(3),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportSolution(example.T, std::to_string(3),
                                     folder,
                                     "T");
        // Tikhonov
        example.tikhonov("GCV", example.sigma);
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
        example.parameterizedBC_postProcess(folder, innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    if (0)
    {
        Info << endl;
        Info << "*********************************************************" << endl;
        Info << "Testing the implementation of regularizing techniques for parametrized BC inverse solver"
             <<
             endl;
        word folder = "./ITHACAoutput/test/";
        List<List<scalar>> heatFluxWeights;
        Eigen::VectorXd residualNorms;
        scalar innerField = 1.0;
        double err_lev = 2; //2% measurements error
        example.solveTrue();
        ITHACAstream::exportSolution(example.T, std::to_string(1),
                                     folder,
                                     "Ttrue");
        example.set_gParametrized("rbf", 0.7);
        example.parameterizedBCoffline();
        example.measError(err_lev);
        Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
        Eigen::MatrixXd b = example.Theta.transpose() * (example.Tmeas_err +
                            example.addSol);
        example.getLinearSystem(A, b);
        Info << "example.gWeights.size() = " << example.gWeights.size() << endl;
        heatFluxWeights.resize(example.gWeights.size() - 1);
        residualNorms.resize(example.gWeights.size() - 1);

        for (label singValueI = 0; singValueI < example.gWeights.size() - 1;
                singValueI++)
        {
            example.TSVD("None", example.sigma, singValueI + 1);
            volScalarField gParametrizedField = example.list2Field(example.g);
            ITHACAstream::exportSolution(gParametrizedField, std::to_string(singValueI + 1),
                                         folder,
                                         "gParametrized");
            ITHACAstream::exportSolution(example.T, std::to_string(singValueI + 1),
                                         folder,
                                         "T");
            heatFluxWeights[singValueI] = example.gWeights;
            Eigen::MatrixXd res = example.linSys.A * example.linSys.x - example.linSys.b;
            residualNorms(singValueI) = Foam::sqrt((res.transpose() * res)(0, 0));
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(example.linSys.A,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd singularValues = svd.singularValues();
        ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen",
                                   folder);
        ITHACAstream::exportVector(residualNorms, "residuals2norm", "eigen",
                                   folder);
        std::cout << "Residual 2-norm = " << residualNorms << std::endl;
        example.parameterizedBC_postProcess(folder, innerField);
        Info << "*********************************************************" << endl;
        Info << endl;
    }

    return 0;
}



