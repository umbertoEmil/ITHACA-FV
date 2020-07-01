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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedInverseLaplacian class

#include "reducedInverseLaplacian.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedInverseLaplacian::reducedInverseLaplacian()
{
}

reducedInverseLaplacian::reducedInverseLaplacian(inverseLaplacianProblem&
        problem)
    :
    problem(&problem)
{
}

void reducedInverseLaplacian::setSensitivityBCandIF()
{
    fvMesh& mesh = problem->_mesh();
    volScalarField& deltaT = problem->_deltaT();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(deltaT, patchI, problem->homogeneousBCcoldSide,
                                           problem->refGrad,
                                           problem->valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(deltaT, patchI, problem->homogeneousBC);
        }
    }
    ITHACAutilities::assignIF(deltaT, problem->homogeneousBC);
}

void reducedInverseLaplacian::assembleLinSys(volScalarField& field)
{
    List<Eigen::MatrixXd> LinSys;
    fvMesh& mesh = problem->_mesh();
    List<scalar> cellVol = mesh.V();

    if (field.name() == "T")
    {
        problem->g = 0.0;
        problem->assignDirectBC();
        ITHACAutilities::assignIF(field, problem->homogeneousBC);
        autoPtr<fvScalarMatrix> pfieldEqn
        (
            new fvScalarMatrix(fvm::laplacian(problem->DT, field))
        );
        LinSysDirect = problem->Tmodes.project(pfieldEqn(), problem->NmodesT);
        pfieldEqn.clear();
        Preduced = Eigen::MatrixXd::Zero(LinSysDirect[1].rows(), 1);
        PreducedOld = Eigen::MatrixXd::Zero(LinSysDirect[1].rows(), 1);
        Treduced.resize(problem->NmodesT, 1);
    }
    else if (field.name() == "lambda")
    {
        problem->assignAdjointBC();
        ITHACAutilities::assignIF(field, problem->homogeneousBC);
        autoPtr<fvScalarMatrix> pfieldEqn
        (
            new fvScalarMatrix(fvm::laplacian(problem->DT, field))
        );
        LinSysAdjoint = problem->lambdaModes.project(pfieldEqn(),
                        problem->NmodesLambda);
        pfieldEqn.clear();
        lambdaReduced.resize(problem->NmodesLambda, 1);
    }
    else if (field.name() == "deltaT")
    {
        setSensitivityBCandIF();
        autoPtr<fvScalarMatrix> pfieldEqn
        (
            new fvScalarMatrix(fvm::laplacian(problem->DT, field))
        );
        LinSysSensitivity = problem->deltaTmodes.project(pfieldEqn(),
                            problem->NmodesDeltaT);
        pfieldEqn.clear();
        deltaTreduced.resize(problem->NmodesDeltaT, 1);
    }
    else
    {
        Info << "Entered wrong field, " << field.name() << endl;
        Info << "It should be T, lambda or deltaT" << endl;
        Info << "Exiting" << endl;
        exit(10);
    }
}

void reducedInverseLaplacian::solveOnline(volScalarField& field, label halfFull)
{
    List<Eigen::MatrixXd> LinSys;
    Eigen::VectorXd residual;
    fvMesh& mesh = problem->_mesh();
    List<scalar> cellVol = mesh.V();

    if (field.name() == "T")
    {
        if (halfFull)
        {
            problem->assignDirectBC();
            autoPtr<fvScalarMatrix> pfieldEqn
            (
                new fvScalarMatrix(fvm::laplacian(problem->DT, field))
            );
            LinSys = problem->Tmodes.project(pfieldEqn(), problem->NmodesT);
            Treduced.resize(LinSys[1].rows(), 1);
            pfieldEqn.clear();
        }
        else
        {
            if (problem->cgIter > 0)
            {
                PreducedOld = gamma * Preduced;
                Preduced = - Poffline * lambdaReduced + PreducedOld;//gamma * PreducedOld;
                LinSysDirect[1] -= beta / problem->k * Preduced;
            }

            LinSys = LinSysDirect;
        }

        Treduced = solveLinearSys(LinSys, Treduced, residual);
    }
    else if (field.name() == "lambda")
    {
        if (halfFull)
        {
            autoPtr<volScalarField> sourcePtr
            (
                new volScalarField(problem->assignAdjointBCandSource())
            );
            autoPtr<fvScalarMatrix> pfieldEqn
            (
                new fvScalarMatrix(fvm::laplacian(problem->DT, field) == - sourcePtr())
            );
            LinSys = problem->lambdaModes.project(pfieldEqn(), problem->NmodesLambda);
            sourcePtr.clear();
            pfieldEqn.clear();
            lambdaReduced.resize(LinSys[1].rows(), 1);
        }
        else
        {
            LinSys = LinSysAdjoint;
            LinSys[1] = LinSysAdjoint[1] - (F_T * Treduced - TmeasReduced) * problem->k;
        }

        lambdaReduced = solveLinearSys(LinSys, lambdaReduced, residual);
        lambdaReducedList.resize(problem->cgIter + 1);
        lambdaReducedList[problem->cgIter] = lambdaReduced;
    }
    else if (field.name() == "deltaT")
    {
        if (halfFull)
        {
            problem->assignSensitivityBC();
            autoPtr<fvScalarMatrix> pfieldEqn
            (
                new fvScalarMatrix(fvm::laplacian(problem->DT, field))
            );
            LinSys = problem->deltaTmodes.project(pfieldEqn(), problem->NmodesDeltaT);
            pfieldEqn.clear();
            deltaTreduced.resize(LinSys[1].rows(), 1);
        }
        else
        {
            if (problem->cgIter == 0)
            {
                b_deltaT_g = - B_deltaT * lambdaReduced / problem->k;
                LinSysSensitivity[1] += b_deltaT_g;
            }
            else
            {
                b_deltaT_g = - B_deltaT * lambdaReduced / problem->k + gamma * b_deltaT_g;
                LinSysSensitivity[1] = - B_deltaT * lambdaReduced / problem->k + gamma *
                                       LinSysSensitivity[1];
            }

            LinSys = LinSysSensitivity;
        }

        deltaTreduced = solveLinearSys(LinSys, deltaTreduced, residual);
    }
    else
    {
        Info << "Entered wrong field, " << field.name() << endl;
        Info << "It should be T, lambda or deltaT" << endl;
        Info << "Exiting" << endl;
        exit(10);
    }
}

void reducedInverseLaplacian::get_L_T()
{
    if (problem->Tmodes.EigenModes.size() == 0)
    {
        problem->Tmodes.toEigen();
    }

    M_Assert(problem->NmodesT <= problem->Tmodes.EigenModes[0].cols(),
             "Number of required modes for direct problem projection is higher then the number of available ones");

    if (problem->NmodesT == 0)
    {
        L_T  = problem->Tmodes.EigenModes[0];
    }
    else
    {
        L_T  = problem->Tmodes.EigenModes[0].leftCols(problem->NmodesT);
    }
}

void reducedInverseLaplacian::get_L_TinterpolationPlane()
{
    L_TinterpolationPlane = L_T;
    Tmeas.resize(L_T.rows(), 1);
    int cellInterpPlaneID = 0;

    for (int i = 0; i < L_TinterpolationPlane.rows(); i++)
    {
        if (i == problem->interpolationPlane.cellID[cellInterpPlaneID])
        {
            Tmeas(i) = problem->interpolationPlane.Tmeas[cellInterpPlaneID];
            cellInterpPlaneID++;
        }
        else
        {
            Tmeas(i) = 0;
            L_TinterpolationPlane.row(i).setZero();
        }
    }
}

void reducedInverseLaplacian::get_L_lambda()
{
    if (problem->lambdaModes.EigenModes.size() == 0)
    {
        problem->lambdaModes.toEigen();
    }

    M_Assert(problem->NmodesLambda <= problem->lambdaModes.EigenModes[0].cols(),
             "Number of required modes for adjoint problem projection is higher then the number of available ones");

    if (problem->NmodesLambda == 0)
    {
        L_lambda  = problem->lambdaModes.EigenModes[0];
    }
    else
    {
        L_lambda = problem->lambdaModes.EigenModes[0].leftCols(problem->NmodesLambda);
    }
}

void reducedInverseLaplacian::get_L_deltaT()
{
    Modes<scalar>& deltaTmodes = problem->deltaTmodes;

    if (deltaTmodes.EigenModes.size() == 0)
    {
        deltaTmodes.toEigen();
    }

    M_Assert(problem->NmodesDeltaT <= deltaTmodes.EigenModes[0].cols(),
             "Number of required modes for direct problem projection is higher then the number of available ones");

    if (problem->NmodesDeltaT == 0)
    {
        L_deltaT  = deltaTmodes.EigenModes[0];
    }
    else
    {
        L_deltaT  = deltaTmodes.EigenModes[0].leftCols(problem->NmodesDeltaT);
    }
}

void reducedInverseLaplacian::assembleG()
{
    fvMesh& mesh = problem->_mesh();
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(L_lambda.rows(), L_lambda.cols());
    label patchID = problem->hotSide_ind;
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];

        for (label mode = 0; mode < problem->NmodesLambda; mode++)
        {
            temp(faceOwner, mode) += L_lambda_GammaIn(faceI, mode) * faceArea;
        }
    }
    Poffline = L_T.transpose() * temp;
}

void reducedInverseLaplacian::assembleAdjointSourceReducedMatrices()
{
    /// Multiply each element on the rows of L_TinterpolationPlane by the correspondent cell volume
    Eigen::MatrixXd tempT = L_TinterpolationPlane;
    Eigen::MatrixXd weigthedTmeas = Tmeas;
    fvMesh& mesh = problem->_mesh();
    List<scalar> cellVol = mesh.V();

    for (int cellI = 0; cellI < L_TinterpolationPlane.rows(); cellI++)
    {
        for (int mode = 0; mode < L_TinterpolationPlane.cols(); mode++)
        {
            tempT(cellI, mode) *= cellVol[cellI];
        }

        weigthedTmeas(cellI) *= cellVol[cellI];
    }

    F_T = L_lambda.transpose() * tempT;
    TmeasReduced = L_lambda.transpose() * weigthedTmeas;
}

void reducedInverseLaplacian::sensitivityBoundaryCoeff()
{
    fvMesh& mesh = problem->_mesh();
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(L_lambda.rows(), L_lambda.cols());
    label patchID = problem->hotSide_ind;
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    L_lambda_GammaIn.resize(cPatch.size(), L_lambda.cols());
    faceArea_GammaIn = Eigen::MatrixXd::Zero(cPatch.size(), 1);
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        faceArea_GammaIn(faceI) = mesh.magSf().boundaryField()[patchID][faceI];

        for (label mode = 0; mode < problem->NmodesLambda; mode++)
        {
            L_lambda_GammaIn(faceI, mode) =
                problem->lambdaModes.EigenModes[patchID + 1](faceI, mode);
            temp(faceOwner, mode) += problem->lambdaModes.EigenModes[patchID + 1](faceI,
                                     mode)
                                     * faceArea_GammaIn(faceI);
        }
    }
    B_deltaT = L_deltaT.transpose() * temp;
}

void reducedInverseLaplacian::conjugateCoeffOffline()
{
    fvMesh& mesh = problem->_mesh();
    label patchID = problem->hotSide_ind;
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    label NmodesLambda = problem->NmodesLambda;
    Eigen::VectorXd faceAreas;
    faceAreas.resize(cPatch.size());
    Eigen::MatrixXd sqrL_Gamma_in = Eigen::MatrixXd::Zero(cPatch.size(),
                                    NmodesLambda * NmodesLambda);;
    forAll(cPatch, faceI)
    {
        label p = 0;
        faceAreas(faceI) = mesh.magSf().boundaryField()[patchID][faceI];

        for (label modeI = 0; modeI < NmodesLambda; modeI++)
        {
            scalar baseFaceValueI = problem->lambdaModes.EigenModes[patchID + 1](faceI,
                                    modeI);

            for (label modeJ = 0; modeJ < NmodesLambda; modeJ++)
            {
                scalar baseFaceValueJ = problem->lambdaModes.EigenModes[patchID + 1](faceI,
                                        modeJ);
                sqrL_Gamma_in(faceI, p) = baseFaceValueI * baseFaceValueJ;
                p++;
            }
        }
    }
    gammaOffline = sqrL_Gamma_in.transpose() * faceAreas;
}

void reducedInverseLaplacian::searchStepOffline()
{
    betaNum1 = Eigen::MatrixXd::Zero(1, problem->NmodesT * problem->NmodesDeltaT);
    betaDen = Eigen::MatrixXd::Zero(problem->NmodesDeltaT * problem->NmodesDeltaT,
                                    1);
    label numberOfCells = L_T.rows();
    fvMesh& mesh = problem->_mesh();
    List<scalar> cellVol = mesh.V();
    label cellInPlaneID = 0;
    Eigen::MatrixXd TmeasWeigthed = Tmeas;

    for (label cellI = 0; cellI < numberOfCells; cellI++)
    {
        label i = 0;
        label j = 0;

        if (cellI == problem->interpolationPlane.cellID[cellInPlaneID])
        {
            TmeasWeigthed(cellI) *= cellVol[cellI];
            cellInPlaneID++;

            for (label col = 0; col < problem->NmodesT * problem->NmodesDeltaT; col++)
            {
                if (j == problem->NmodesDeltaT)
                {
                    j = 0;
                    i++;
                }

                betaNum1(0, col) += cellVol[cellI] * L_T(cellI, i) * L_deltaT(cellI, j);
                j++;
            }

            label p = 0;

            for (label modeI = 0; modeI < problem->NmodesDeltaT; modeI++)
            {
                for (label modeJ = 0; modeJ < problem->NmodesDeltaT; modeJ++)
                {
                    betaDen(p) += L_deltaT(cellI, modeI) * L_deltaT(cellI, modeJ) * cellVol[cellI];
                    p++;
                }
            }
        }
    }

    betaNum2 = TmeasWeigthed.transpose() * L_deltaT;
}

void reducedInverseLaplacian::computeConjugateCoeff()
{
    Eigen::VectorXd lambdaReducedSqrOld = lambdaReducedSqr;

    if (problem->cgIter == 0)
    {
        lambdaReducedSqr.resize(problem->NmodesLambda * problem->NmodesLambda);
    }

    for (label i = 0; i < problem->NmodesLambda; i++)
    {
        for (label j = 0; j < problem->NmodesLambda; j++)
        {
            lambdaReducedSqr(i * problem->NmodesLambda + j) = lambdaReduced(
                        i) * lambdaReduced(j);
        }
    }

    if (problem->cgIter == 0)
    {
        gamma = lambdaReducedSqr.dot(gammaOffline);
    }
    else
    {
        gamma = lambdaReducedSqr.dot(gammaOffline) / lambdaReducedSqrOld.dot(
                    gammaOffline);
    }

    gammaList.resize(problem->cgIter + 1);
    gammaList[problem->cgIter] = gamma;
    std::cout << "gamma = " << gamma << std::endl;
}

void reducedInverseLaplacian::searchStep()
{
    Eigen::MatrixXd v_TdeltaT;
    Eigen::MatrixXd tempNum;
    Eigen::MatrixXd tempDen;
    Eigen::VectorXd sqrDeltaTreduced;
    v_TdeltaT.resize(problem->NmodesT * problem->NmodesDeltaT, 1);
    sqrDeltaTreduced.resize(problem->NmodesDeltaT * problem->NmodesDeltaT);

    for (label i = 0; i < problem->NmodesT; i++)
    {
        for (label j = 0; j < problem->NmodesDeltaT; j++)
        {
            v_TdeltaT(j + i * problem->NmodesDeltaT, 0) = Treduced(i, 0) * deltaTreduced(j,
                    0);
        }
    }

    for (label i = 0; i < problem->NmodesDeltaT; i++)
    {
        for (label j = 0; j < problem->NmodesDeltaT; j++)
        {
            sqrDeltaTreduced(i * problem->NmodesDeltaT + j) = deltaTreduced(
                        i) * deltaTreduced(j);
        }
    }

    tempNum = betaNum1 * v_TdeltaT - betaNum2 * deltaTreduced;
    beta = tempNum(0, 0) / betaDen.dot(sqrDeltaTreduced);
    betaList.resize(problem->cgIter + 1);
    betaList[problem->cgIter] = beta;
    Info << "beta = " << beta << endl;
}

void reducedInverseLaplacian::reconstructHeatFlux()
{
    g = Eigen::MatrixXd::Zero(problem->g.size(), 1);
    P.conservativeResize(g.rows());
    forAll (lambdaReducedList, iter)
    {
        P = - L_lambda_GammaIn * lambdaReducedList[iter] + gammaList[iter] * P;
        g -= betaList[iter] * P;
    }
}

void reducedInverseLaplacian::writeFields(label folderNumber,
        const char* folder)
{
    fvMesh& mesh = problem->_mesh();
    volScalarField& T = problem->_T();
    volScalarField& lambda = problem->_lambda();
    volScalarField& deltaT = problem->_deltaT();
    Info << "Code well the fields reconstruction\n";
    exit(10);
    //T = problem->Tmodes.reconstruct(Treduced, "T");
    //lambda = problem->lambdaModes.reconstruct(lambdaReduced, "lambda");
    //deltaT = problem->deltaTmodes.reconstruct(deltaTreduced, "deltaT");
    reconstructHeatFlux();
    autoPtr<volScalarField> gVolField_
    (
        new volScalarField("g", T)
    );
    volScalarField& gVolField = gVolField_();
    ITHACAutilities::assignIF(gVolField, problem->homogeneousBC);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[problem->hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        gVolField[faceOwner] = g(faceI);
    }
    ITHACAstream::exportSolution(T, std::to_string(folderNumber + 1), folder,
                                 "T");
    ITHACAstream::exportSolution(lambda, std::to_string(folderNumber + 1), folder,
                                 "lambda");
    ITHACAstream::exportSolution(deltaT, std::to_string(folderNumber + 1), folder,
                                 "deltaT");
    ITHACAstream::exportSolution(gVolField, std::to_string(folderNumber + 1),
                                 folder, "g");
}

int reducedInverseLaplacian::conjugateGradientConvergenceCheck()
{
    fvMesh& mesh = problem->_mesh();
    double Jold = J;
    J = 0;
    Eigen::VectorXd sqTdiff = L_TinterpolationPlane * Treduced - Tmeas;
    sqTdiff = sqTdiff.cwiseProduct(sqTdiff);
    forAll(mesh.V(),
           cellI)//this can become a loop only over the intepolation plane cells
    {
        J += 0.5 * sqTdiff(cellI) * mesh.V()[cellI];
    }
    Info << "J =" << J << endl;

    if (J <= problem->Jtol)
    {
        Info << "Convergence reached in " << problem->cgIter << " iterations" << endl;
        return (1);
    }
    else if (Foam::mag((Jold - J) / J) <= problem->JtolRel)
    {
        Info << "Relative tolerance criteria meet in " << problem->cgIter <<
             " iterations" <<
             endl;
        Info << "|Jold - J| / |J| = " << Foam::mag((Jold - J) / J) << endl;
        return (1);
    }
    else if (Foam::mag((Jold - J) / J) > 1e2)
    {
        Info << "J =" << J << endl;
        Info << "The algorithm is diverging. Exiting" << endl;
        return (2);
    }
    else
    {
        return (0);
    }
}

void reducedInverseLaplacian::conjugateGradientOffline()
{
    if (!offlineCG)
    {
        Info << "Performing offline computation of the matrices" << endl;
        volScalarField& T = problem->_T();
        volScalarField& lambda = problem->_lambda();
        volScalarField& deltaT = problem->_deltaT();
        problem->set_g();
        problem->set_valueFraction();
        assembleLinSys(T);
        assembleLinSys(lambda);
        assembleLinSys(deltaT);
        get_L_T();
        get_L_TinterpolationPlane();
        get_L_lambda();
        get_L_deltaT();
        assembleAdjointSourceReducedMatrices();
        sensitivityBoundaryCoeff();
        conjugateCoeffOffline();
        searchStepOffline();
        assembleG();
        offlineCG = 1;
    }
}

int reducedInverseLaplacian::conjugateGradient()
{
    volScalarField& T = problem->_T();
    volScalarField& lambda = problem->_lambda();
    volScalarField& deltaT = problem->_deltaT();
    conjugateGradientOffline();
    problem->cgIter = 0;
    gamma = 0.0;
    assembleLinSys(T);
    Treduced = Eigen::MatrixXd::Zero(problem->NmodesT, 1);
    lambdaReduced = Eigen::MatrixXd::Zero(problem->NmodesLambda, 1);
    deltaTreduced = Eigen::MatrixXd::Zero(problem->NmodesDeltaT, 1);

    /*
    /// RBF in not working yet because I do not know how
    /// to compute the matrix containing the basis of the RBF
    Eigen::MatrixXd& RBFweightsP(RBFweights);
    Eigen::MatrixXd& RBFbasisP(RBFbasis);
    problem->thermocouplesInterpolation(RBFweightsP, RBFbasisP);
    */

    while (problem->cgIter < problem->cgIterMax)
    {
        Info << "Iteration " << problem->cgIter + 1 << endl;
        solveOnline(T);
        label CGoutput = conjugateGradientConvergenceCheck();

        if (CGoutput == 1)
        {
            return (1);
        }
        else if (CGoutput)
        {
            return (2);
        }

        solveOnline(lambda);
        computeConjugateCoeff();
        solveOnline(deltaT);
        searchStep();
        problem->cgIter++;
    }

    return (0);
}

int reducedInverseLaplacian::conjugateGradientHalfFull()
{
    label halfFull = 1;
    volScalarField& T = problem->_T();
    volScalarField& lambda = problem->_lambda();
    volScalarField& deltaT = problem->_deltaT();
    problem->set_g();
    problem->set_valueFraction();
    problem->cgIter = 0;
    problem->J = 0;
    problem->P = problem->g;
    problem->gradJ = problem->g;
    problem->gamma = 0.0;
    problem->gamma_den = 0.0;

    while (problem->cgIter < problem->cgIterMax)
    {
        Info << "Iteration " << problem->cgIter + 1 << endl;
        solveOnline(T, halfFull);
        Info << "Code well the fields reconstruction\n";
        exit(10);
        //T = problem->Tmodes.reconstruct(Treduced, "T");
        problem->differenceBetweenDirectAndMeasure();

        if (problem->conjugateGradientConvergenceCheck())
        {
            return (1);
        }

        solveOnline(lambda, halfFull);
        //lambda = problem->lambdaModes.reconstruct(lambdaReduced, "lambda");
        problem->computeGradJ();
        problem->searchDirection();
        solveOnline(deltaT, halfFull);
        //deltaT = problem->deltaTmodes.reconstruct(deltaTreduced, "deltaT");
        problem->sensibilitySolAtThermocouplesLocations();
        problem->computeSearchStep();
        problem->updateHeatFlux();
        problem->cgIter += 1;
        problem->writeFields(problem->cgIter, "./ITHACAoutput/online_rec");
    }

    return (0);
}

void reducedInverseLaplacian::get_heatFluxBases()
{
    // I have to rewrite the bases for the cells and not for theboundary faces
    Eigen::MatrixXd temp =
        ITHACAstream::readMatrix("./ITHACAoutput/podMarquardt/gReducedBases_mat.txt");
    problem->gWeights.resize(temp.cols());
    fvMesh& mesh = problem->_mesh();
    heatFluxBases = Eigen::MatrixXd::Zero(L_T.rows(), temp.cols());
    label patchID = problem->hotSide_ind;
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        scalar faceArea = mesh.magSf().boundaryField()[patchID][faceI];

        for (label mode = 0; mode < temp.cols(); mode++)
        {
            heatFluxBases(faceOwner, mode) += temp(faceI, mode)
                                              * faceArea;// / problem->k;
        }
    }
}


int reducedInverseLaplacian::Marquardt()
{
    volScalarField& T = problem->_T();
    label iter = 0;
    problem->set_g();
    problem->set_valueFraction();
    assembleLinSys(T);
    Eigen::VectorXd dummy;
    get_L_T();
    Info << "L_T shape = " << L_T.rows() << ", " << L_T.cols() << endl;
    get_heatFluxBases();
    Info << "heatFluxBases shape = " << heatFluxBases.rows() << ", " <<
         heatFluxBases.cols() << endl;
    Eigen::MatrixXd heatFluxBasesReduced = L_T.transpose() * heatFluxBases;
    problem->set_gParametrized("pod");
    problem->MarquardtMethodSetUp();
    // Marquardt parameters
    label iterMax = 100;
    scalar tol = 1e-3;
    // Damping factor update parameters
    scalar dampingFactor = 0; // 1e-5;
    scalar c = 1.1;
    scalar d = 2;
    scalar csi = 2;
    Eigen::VectorXd residual;
    residual.resize(problem->Tmeas.size());
    Eigen::VectorXd Tcomp;
    Tcomp.resize(problem->Tmeas.size());
    Eigen::VectorXd TcompOld = Tcomp;
    Eigen::VectorXd TcompNew = Tcomp;
    List<scalar> gWeightsOld = problem->gWeights;
    List<scalar> gWeightsNew = problem->gWeights;
    Eigen::VectorXd weigthsUpdate;
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    Eigen::MatrixXd Jvector;
    Eigen::MatrixXd JacobianOld = problem->Jacobian;
    label convergence = 0;
    Eigen::VectorXd gWeightsEigen = Eigen::MatrixXd::Zero(problem->gWeights.size(),
                                    1);
    List<Eigen::MatrixXd> LinSysDirectOrig = LinSysDirect;
    word folder = "./ITHACAoutput/testRedMarquardt/";
    auto t1 = std::chrono::high_resolution_clock::now();

    while (iter < iterMax && convergence == 0)
    {
        Info << endl << endl ;
        Info << "Iteration = " << iter + 1 << endl;
        Info << "damping factor = " << dampingFactor << endl;
        std::cout << "gWeights = " << gWeightsEigen << std::endl;
        problem->update_gParametrized(problem->gWeights);
        volScalarField gParametrizedField = problem->list2Field(problem->g);
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(iter + 1),
                                     folder,
                                     "gParametrized");
        // Solve reduced direct problem
        LinSysDirect[1] =  LinSysDirectOrig[1] + heatFluxBasesReduced * gWeightsEigen /
                           problem->k;
        Treduced = solveLinearSys(LinSysDirect, Treduced, residual);
        Info << "Code well the fields reconstruction\n";
        exit(10);
        //T = problem->Tmodes.reconstruct(Treduced, "T");
        ITHACAstream::exportSolution(T, std::to_string(iter + 1),
                                     "./ITHACAoutput/testRedMarquardt/", "T");
        forAll(problem->thermocouplesCellID, cellI)
        {
            Tcomp(cellI) = (L_T * Treduced)(problem->thermocouplesCellID[cellI], 0);
        }
        // Uncomment for solving full direct
        //problem->solveDirect();
        //forAll(problem->thermocouplesCellID, cellI)
        //{
        //    Tcomp(cellI) = T.internalField()[problem->thermocouplesCellID[cellI]];
        //}
        //ITHACAstream::exportSolution(T, std::to_string(iter + 1), folder,
        //                             "Tdirect");
        residual = Tcomp - problem->Tmeas;
        problem->J = 0.5 * residual.dot(residual);
        Eigen::VectorXd resNonDim = residual.array() / problem->Tmeas.array();
        scalar JnonDim = Foam::sqrt((resNonDim.dot(resNonDim)) / resNonDim.size());
        Info << "J = " << problem->J << endl;
        Info << "JnonDim = " << JnonDim << endl;
        dampingFactor /= d;
        linSys[0] = problem->Jacobian.transpose() * problem->Jacobian + dampingFactor *
                    problem->diagJacobian; //Eigen::MatrixXd::Identity(Jacobian.cols(), Jacobian.cols());
        linSys[1] = - problem->Jacobian.transpose() * residual;
        weigthsUpdate = linSys[0].fullPivLu().solve(linSys[1]);
        forAll(problem->gWeights, weightI)
        {
            problem->gWeights[weightI] += weigthsUpdate(weightI);
        }
        gWeightsEigen += weigthsUpdate;
        scalar gWeightsNorm = gWeightsEigen.transpose() * gWeightsEigen;
        gWeightsNorm = Foam::sqrt(gWeightsNorm);
        Info << "weigthsUpdate.norm() / gWeightsNorm = " << weigthsUpdate.norm() /
             gWeightsNorm << endl;

        if (JnonDim < 0.02
                ||  weigthsUpdate.norm() / gWeightsNorm < tol) // Convergence check
        {
            Info << "Levenberg-Marquardt algorithm converged" << endl;
            convergence = 1;
        }

        iter++;
        //        TcompNew = Tcomp;
        //        residual = Tcomp - Tmeas;
        //        //std::cout << "residual = " << residual << std::endl;
        //        scalar Jold = J;
        //        J = 0.5 * residual.dot(residual);
        //        Eigen::VectorXd resNonDim = residual.array() / Tmeas.array();
        //        scalar JnonDim = Foam::sqrt((resNonDim.dot(resNonDim)) / resNonDim.size());
        //        Info << "J = " << J << endl;
        //        Info << "JnonDim = " << JnonDim << endl;
        //
        //        if (exportSolutions)
        //        {
        //            Jvector.conservativeResize(iter + 1, 1);
        //            Jvector(iter) = J;
        //        }
        //
        //        //dampingFactor = JnonDim * JnonDim * JnonDim * JnonDim;
        //        dampingFactor /= d;
        //
        //        if (iter > 0 && updateJacobian)
        //        {
        //            Info << "Updating Jacobian" << endl;
        //            Jacobian = MarquardtUpdateJacobian(gWeightsOld, weigthsUpdate, Tcomp, TcompOld);
        //        }
        //
        //        linSys[0] = Jacobian.transpose() * Jacobian + dampingFactor *
        //                    diagJacobian; //Eigen::MatrixXd::Identity(Jacobian.cols(), Jacobian.cols());
        //        linSys[1] = - Jacobian.transpose() * residual;
        //        weigthsUpdate = linSys[0].fullPivLu().solve(linSys[1]);
        //        gWeightsOld = gWeights;
        //        TcompOld = TcompNew;
        //        forAll(gWeights, weigthI)
        //        {
        //            gWeights[weigthI] += weigthsUpdate(weigthI);
        //        }
        //        scalar gWeightsNorm = 0;
        //        forAll(gWeights, weigthI)
        //        {
        //            gWeightsNorm += gWeights[weigthI] * gWeights[weigthI];
        //        }
        //        gWeightsNorm = Foam::sqrt(gWeightsNorm);
        //        Info << "weigthsUpdate.norm() / gWeightsNorm = " << weigthsUpdate.norm() /
        //             gWeightsNorm << endl;
        //
        //        if (JnonDim < 0.02
        //                ||  weigthsUpdate.norm() / gWeightsNorm < tol) // Convergence check
        //        {
        //            Info << "Levenberg-Marquardt algorithm converged" << endl;
        //            convergence = 1;
        //        }
        //
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>
                     (t2 - t1);
    double time = time_span.count();
    std::cout << "CPU time reduced order = " << time << std::endl;
    return (0);
}

// ************************************************************************* //

