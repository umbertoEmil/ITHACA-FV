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
/// Source file of the inverseLaplacianProblem class.


#include "inverseLaplacianProblem.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
inverseLaplacianProblem::inverseLaplacianProblem() {}

inverseLaplacianProblem::inverseLaplacianProblem(int argc, char* argv[])
    :
    DT("DT", dimensionSet(1, 1, -3, -1, 0, 0, 0), 1.0)
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createThermocouples.H"
    thermocouplesPos = TCpos;
#include "createFvOptions.H"
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    nProcs = Pstream::nProcs();
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void inverseLaplacianProblem::set_g()
{
    volScalarField& T = _T();
    g.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
    }
}

void inverseLaplacianProblem::set_gBaseFunctions(word type,
        scalar shapeParameter)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();

    if (type == "rbf")
    {
        Info << "Radial Basis Functions are used." << endl;
        // The center of each function is the projection of each thermocouple
        // on the boundary hotSide

        if (!thermocouplesRead)
        {
            readThermocouples();
        }

        gBaseFunctions.resize(thermocouplesNum);
        gWeights.resize(thermocouplesNum);
        //for (label funcI = 0; funcI < thermocouplesNum; funcI++)
        forAll(gBaseFunctions, funcI)
        {
            gBaseFunctions[funcI].resize(T.boundaryField()[hotSide_ind].size());
            scalar thermocoupleX = thermocouplesPos[funcI].x();
            scalar thermocoupleZ = thermocouplesPos[funcI].z();
            forAll (T.boundaryField()[hotSide_ind], faceI)
            {
                scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - thermocoupleX) +
                                           (faceZ - thermocoupleZ) * (faceZ - thermocoupleZ));
                gBaseFunctions[funcI][faceI] = Foam::exp(- (shapeParameter *
                                               shapeParameter
                                               * radius * radius));
            }
        }
    }
    else if (type == "pod")
    {
        Eigen::MatrixXd temp =
            ITHACAstream::readMatrix("./ITHACAoutput/podMarquardt/gReducedBases_mat.txt");
        gBaseFunctions.resize(temp.cols());
        gWeights.resize(temp.cols());
        forAll(gBaseFunctions, baseI)
        {
            gBaseFunctions[baseI].resize(temp.rows());
            forAll(gBaseFunctions[baseI], faceI)
            {
                gBaseFunctions[baseI][faceI] = temp(faceI, baseI);
            }
        }
    }
}

void inverseLaplacianProblem::set_gBaseFunctionsPOD(label Nmodes)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    Eigen::MatrixXd gBaseFuncEigen;
    set_gParametrized("rbf");

    if (Nmodes == 0)
    {
        Info << "Selecting all modes." << endl;
        Nmodes = gBaseFunctions.size();
    }

    gBaseFuncEigen.resize(gBaseFunctions[0].size(), gBaseFunctions.size());
    Eigen::VectorXd faceAreaVect;
    faceAreaVect.resize(mesh.magSf().boundaryField()[hotSide_ind].size());
    forAll(gBaseFunctions, funcI)
    {
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            if (funcI == 0)
            {
                faceAreaVect(faceI) = mesh.magSf().boundaryField()[hotSide_ind][faceI];
            }

            gBaseFuncEigen(faceI, funcI) = gBaseFunctions[funcI][faceI];
        }
    }
    Eigen::MatrixXd correlationMatrix = gBaseFuncEigen.transpose() *
                                        faceAreaVect.asDiagonal() * gBaseFuncEigen;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(correlationMatrix,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    gPODmodes = svd.matrixU().leftCols(Nmodes);
    Eigen::MatrixXd gBaseFuncEigen_new = gBaseFuncEigen * gPODmodes;
    Info << "gBaseFuncEigen_new size = " << gBaseFuncEigen_new.cols() << ", " <<
         gBaseFuncEigen_new.rows() << endl;
    gBaseFunctions.resize(Nmodes);
    gWeights.resize(Nmodes);
    forAll(gBaseFunctions, funcI)
    {
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            gBaseFunctions[funcI][faceI] = gBaseFuncEigen_new(faceI, funcI);
        }
    }
    //// Resizing Theta
    //if(Theta.cols() > 0)
    //{
    //    Theta = Theta.leftCols(Nmodes);
    //}
    //else
    //{
    //    Info << "Theta not computed yet" << endl;
    //}
}

void inverseLaplacianProblem::set_gParametrized(word baseFuncType,
        scalar shapeParameter)
{
    volScalarField& T = _T();
    set_gBaseFunctions(baseFuncType, shapeParameter);
    g.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    forAll (gWeights, weigthI)
    {
        gWeights[weigthI] = 0; //-10000;
    }
    Info << "gWeights = " << gWeights << endl;
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
        forAll (gWeights, weigthI)
        {
            g[faceI] += gWeights[weigthI] * gBaseFunctions[weigthI][faceI];
        }
    }
}

void inverseLaplacianProblem::update_gParametrized(List<scalar> weigths)
{
    M_Assert(weigths.size() == gBaseFunctions.size(),
             "weigths size different from basis functions size");
    volScalarField& T = _T();
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        g[faceI] = 0.0;
        forAll (weigths, weigthI)
        {
            g[faceI] += weigths[weigthI] * gBaseFunctions[weigthI][faceI];
        }
    }
}

volScalarField inverseLaplacianProblem::list2Field(List<scalar> list,
        scalar innerField)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    volScalarField field(T);
    ITHACAutilities::assignIF(field, innerField);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        field[faceOwner] = list[faceI];
    }
    return field;
}

volScalarField inverseLaplacianProblem::eigen2Field(Eigen::MatrixXd matrix)
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    volScalarField field(T);
    ITHACAutilities::assignIF(field, homogeneousBC);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        field[faceOwner] = matrix(faceI);
    }
    return field;
}

Eigen::MatrixXd inverseLaplacianProblem::MarquardtUpdateJacobian(
    List<scalar> gWeightsOld,
    Eigen::VectorXd weigthsUpdate, Eigen::VectorXd Tdirect,
    Eigen::VectorXd TdirectOld)
{
    volScalarField& T = _T();
    List<scalar> gWeights = gWeightsOld;
    Eigen::MatrixXd Jacobian;
    Jacobian.resize(Tmeas.size(), gWeights.size());

    for (label j = 0; j < Jacobian.cols(); j++)
    {
        gWeights[j] = gWeightsOld[j] + weigthsUpdate(j);
        update_gParametrized(gWeights);
        solveDirect();
        Tdirect = fieldValueAtThermocouples(T);

        for (label i = 0; i < Jacobian.rows(); i++)
        {
            Jacobian(i, j) = (Tdirect(i) - TdirectOld(i)) / (gWeights[j] - gWeightsOld[j]);
        }

        gWeights = gWeightsOld;
    }

    return Jacobian;
}

Eigen::VectorXd  inverseLaplacianProblem::TSVD(Eigen::MatrixXd A,
        Eigen::MatrixXd b, label filter)
{
    // Add check on b
    Info << "Using truncated SVD for regularization" << endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd x;

    for (label i = 0; i < filter; i++)
    {
        scalar coeff = (U.col(i).transpose() * b)(0, 0);

        if (i == 0)
        {
            // / svd.singularValues()(i) << endl;
            x = coeff / svd.singularValues()(i) * V.col(i);
        }
        else
        {
            x += coeff / svd.singularValues()(i) * V.col(i);
        }
    }

    return x;
}

void inverseLaplacianProblem::parameterizedBCoffline(bool force)
{
    fvMesh& mesh = _mesh();
    Tbasis.resize(0);
    Tad_base.resize(0);

    if (ITHACAutilities::check_file(folderOffline + "Theta_mat.txt") && force == 0)
    {
        Info << "\nOffline already computed." << endl;
	Info << "Check that the basis used for the parameterized BC are correct (RBF, POD, etc.)\n";
        Theta = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");
	addSol = ITHACAstream::readMatrix(folderOffline + "addSol_mat.txt");

        volScalarField& T(_T());
        ITHACAstream::read_fields(Tad_base, "Tad", folderOffline, 0, 1);

        ITHACAstream::read_fields(Tbasis, "T",
                                  folderOffline);
    }
    else
    {
        Info << "\nComputing offline" << endl;
        solveAdditional();
	ITHACAstream::exportVector(addSol, "addSol", "eigen", folderOffline);
	M_Assert(Tmeas.size() > 0, "Initialize Tmeas");
	M_Assert(gWeights.size() > 0, "Initialize gWeights");
        Theta.resize(Tmeas.size(), gWeights.size());

        for (label j = 0; j < Theta.cols(); j++)
        {
            gWeights = Foam::zero();
            gWeights[j] =  1;
            update_gParametrized(gWeights);
            Info << "Solving for j = " << j << endl;
            solveDirect();
            volScalarField& T = _T();
	    Tbasis.append(T);
            Tdirect = fieldValueAtThermocouples(T);

            for (label i = 0; i < Theta.rows(); i++)
            {
                Theta(i, j) = Tdirect(i) + addSol(i);
            }

            volScalarField gParametrizedField = list2Field(g);
            ITHACAstream::exportSolution(gParametrizedField, std::to_string(j + 1),
                                         folderOffline,
                                         "gParametrized");
        }
	ITHACAstream::exportFields(Tbasis, folderOffline, "T");

        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
        Info << "\nOffline part ENDED\n" << endl;
    }

    Eigen::MatrixXd A = Theta.transpose() * Theta;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues = svd.singularValues();
    double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
    Info << "Condition number = " << conditionNumber << endl;
    ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen",
                               folderOffline);
}

Eigen::VectorXd inverseLaplacianProblem::parameterizedBC(word linSys_solver,
        label TSVD_filter)
{
    Info << endl << "Using quasilinearity of direct problem" << endl;
    //parameterizedBCoffline(folder, forceOffline);
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    Info << "Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;
    linSys[0] = Theta.transpose() * Theta;
    linSys[1] = Theta.transpose() * (Tmeas + addSol);
    Eigen::VectorXd weigths;

    if (linSys_solver == "fullPivLU")
    {
        weigths = linSys[0].fullPivLu().solve(linSys[1]);
    }
    else if (linSys_solver == "jacobiSvd")
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        weigths = svd.solve(linSys[1]);
    }
    else if (linSys_solver == "householderQr")
    {
        weigths = linSys[0].householderQr().solve(linSys[1]);
    }
    else if (linSys_solver == "ldlt")
    {
        weigths = linSys[0].ldlt().solve(linSys[1]);
    }
    else if (linSys_solver == "inverse")
    {
        weigths = linSys[0].inverse() * linSys[1];
    }
    else if (linSys_solver == "TSVD")
    {
        weigths = TSVD(linSys[0], linSys[1], TSVD_filter);
    }
    else
    {
        Info << "Select a linear system solver in this list:" << endl
             << "fullPivLU, jacobiSvd, householderQr, ldlt, inverse, TSVD" << endl;
        exit(1);
    }
    parameterizedBCpostProcess(weigths);
    return weigths;
}

void inverseLaplacianProblem::parameterizedBCpostProcess(Eigen::VectorXd weigths)
{
    //// Printing outputs at screen
    //std::cout << "Eigenvalues of Theta.transpose() * Theta are " << std::endl;
    //std::cout << linSys[0].eigenvalues() << std::endl;
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0], Eigen::ComputeThinU | Eigen::ComputeThinV);
    //std::cout << "Singular values of Theta.transpose() * Theta are " << std::endl;
    //std::cout << svd.singularValues() << std::endl;
    //std::cout << "debug: weigths size = " << std::endl;
    //std::cout << weigths.size() << std::endl;
    //
    //residual =  linSys[0] * weigths - linSys[1];
    //std::cout << "Residual  = " << std::endl;
    //std::cout << residual << std::endl;
    //std::cout << "Residual 2-norm = " << std::endl;
    //std::cout << residual.squaredNorm() << std::endl;
    gWeights.resize(weigths.size());
    forAll(gWeights, weightI)
    {
        gWeights[weightI] = weigths(weightI);
    }
    update_gParametrized(gWeights);

    reconstructT();
    volScalarField& T = _T();
    Tdirect = fieldValueAtThermocouples(T);

    J = 0.5 * (Tdirect - Tmeas).dot(Tdirect - Tmeas);
    Info << "J = " << J << endl;
    Info << "End" << endl;
    Info << endl;
}

void inverseLaplacianProblem::MarquardtMethodSetUp()
{
    Info << endl << "Setting up the Marquardt method" << endl;
    Jacobian.resize(Tmeas.size(), gWeights.size());
    List<scalar> gWeightsOld = gWeights;
    Eigen::VectorXd TdirectOld = Tdirect;
    //Initializing the Jacobian
    update_gParametrized(gWeights);
    solveDirect();
    volScalarField& T = _T();
    Tdirect = fieldValueAtThermocouples(T);

    for (label j = 0; j < Jacobian.cols(); j++)
    {
        if (Foam::mag(gWeights[j]) > 1e-8)
        {
            gWeightsOld[j] = gWeights[j] + gWeights[j] / 10;
        }
        else
        {
            gWeightsOld[j] = gWeights[j] - 1e5 ;//+ 1e3;//e5;
        }

        update_gParametrized(gWeightsOld);
        Info << "Solving for j = " << j << endl;
        solveDirect();
        Tdirect = fieldValueAtThermocouples(T);

        for (label i = 0; i < Jacobian.rows(); i++)
        {
            Jacobian(i, j) = (Tdirect(i) - TdirectOld(i)) / (gWeights[j] - gWeightsOld[j]);
        }

        gWeightsOld = gWeights;
    }

    std::cout << "Jacobian = " << Jacobian << std::endl;
    Eigen::MatrixXd diag = (Jacobian.transpose() * Jacobian).diagonal();
    diagJacobian = diag.matrix().asDiagonal();
    Info << "End of the set up" << endl;
    Info << endl;
}

void inverseLaplacianProblem::MarquardtMethod(label exportSolutions,
        word folder,
        label updateJacobian)
{
    volScalarField& T = _T();
    // Marquardt parameters
    label iter = 0;
    label iterMax = 100;
    scalar tol = 1e-3;
    // Damping factor update parameters
    scalar dampingFactor = 10;
    scalar c = 1.1;
    scalar d = 2;
    scalar csi = 2;
    Eigen::VectorXd residual;
    residual.resize(Tmeas.size());
    Eigen::VectorXd TdirectOld = Tdirect;
    Eigen::VectorXd TdirectNew = Tdirect;
    List<scalar> gWeightsOld = gWeights;
    List<scalar> gWeightsNew = gWeights;
    Eigen::VectorXd weigthsUpdate;
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    Eigen::MatrixXd Jvector;
    Eigen::MatrixXd JacobianOld = Jacobian;
    label convergence = 0;
    auto t1 = std::chrono::high_resolution_clock::now();

    while (iter < iterMax && convergence == 0)
    {
        Info << endl << endl ;
        Info << "Iteration = " << iter + 1 << endl;
        Info << "damping factor = " << dampingFactor << endl;
        update_gParametrized(gWeights);
        Info << "gWeights = " << gWeights << endl;
        volScalarField gParametrizedField = list2Field(g);
        solveDirect();

        if (exportSolutions)
        {
            ITHACAstream::exportSolution(gParametrizedField, std::to_string(iter + 1),
                                         folder,
                                         "gParametrized");
            ITHACAstream::exportSolution(T, std::to_string(iter + 1), folder,
                                         "T");
        }

        Tdirect = fieldValueAtThermocouples(T);
        TdirectNew = Tdirect;
        residual = Tdirect - Tmeas;
        //std::cout << "residual = " << residual << std::endl;
        scalar Jold = J;
        J = 0.5 * residual.dot(residual);
        Eigen::VectorXd resNonDim = residual.array() / Tmeas.array();
        scalar JnonDim = Foam::sqrt((resNonDim.dot(resNonDim)) / resNonDim.size());
        Info << "J = " << J << endl;
        Info << "JnonDim = " << JnonDim << endl;

        if (exportSolutions)
        {
            Jvector.conservativeResize(iter + 1, 1);
            Jvector(iter) = J;
        }

        //dampingFactor = JnonDim * JnonDim * JnonDim * JnonDim;
        dampingFactor /= d;

        if (iter > 0 && updateJacobian)
        {
            Info << "Updating Jacobian" << endl;
            Jacobian = MarquardtUpdateJacobian(gWeightsOld, weigthsUpdate, Tdirect,
                                               TdirectOld);
        }

        linSys[0] = Jacobian.transpose() * Jacobian + dampingFactor *
                    diagJacobian; //Eigen::MatrixXd::Identity(Jacobian.cols(), Jacobian.cols());
        linSys[1] = - Jacobian.transpose() * residual;
        weigthsUpdate = linSys[0].fullPivLu().solve(linSys[1]);
        gWeightsOld = gWeights;
        TdirectOld = TdirectNew;
        forAll(gWeights, weigthI)
        {
            gWeights[weigthI] += weigthsUpdate(weigthI);
        }
        scalar gWeightsNorm = 0;
        forAll(gWeights, weigthI)
        {
            gWeightsNorm += gWeights[weigthI] * gWeights[weigthI];
        }
        gWeightsNorm = Foam::sqrt(gWeightsNorm);
        Info << "weigthsUpdate.norm() / gWeightsNorm = " << weigthsUpdate.norm() /
             gWeightsNorm << endl;

        if (JnonDim < 0.02
                ||  weigthsUpdate.norm() / gWeightsNorm < 1e-2) // Convergence check
        {
            Info << "Levenberg-Marquardt algorithm converged" << endl;
            convergence = 1;
            Info << "Weights old method = " << endl;
            Info << gWeights << endl;
        }

        iter++;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>
                     (t2 - t1);
    double time = time_span.count();
    std::cout << "CPU time Full order = " << time << std::endl;

    if (convergence == 0)
    {
        Info << "Levenberg-Marquardt algorithm did NOT converge" << endl;
    }

    if (exportSolutions)
    {
        update_gParametrized(gWeights);
        volScalarField gParametrizedField = list2Field(g);
        solveDirect();
        ITHACAstream::exportSolution(T, std::to_string(iter + 1), folder,
                                     "T");
        ITHACAstream::exportSolution(gParametrizedField, std::to_string(iter + 1),
                                     folder,
                                     "gParametrized");
        ITHACAstream::exportMatrix(Jvector, "costFunction", "eigen", folder);
    }
}

void inverseLaplacianProblem::set_valueFraction()
{
    fvMesh& mesh = _mesh();
    valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
    homogeneousBCcoldSide.resize(mesh.boundaryMesh()["coldSide"].size());
    valueFractionAdj.resize(mesh.boundaryMesh()["coldSide"].size());
    Eigen::VectorXd faceCellDist =
        ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
    forAll (valueFraction, faceI)
    {
        scalar faceDist = faceCellDist(faceI);
        valueFraction[faceI] =  1.0 / (1.0 + (k / H / faceDist));
        valueFractionAdj[faceI] =  1 / (1 + (1 / k / H / faceDist));
        homogeneousBCcoldSide[faceI] =  0;
    }
    refGrad = homogeneousBCcoldSide;
}


void inverseLaplacianProblem::solveTrue()
{
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    Foam::Time& runTime = _runTime();
    gTrue.resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    scalar moldHeight = 1.2; //[m]
    scalar moldWidth = 1.97; //[m]
    scalar castingSpeed = 0.0375; //[m/s]
    forAll (T.boundaryField()[hotSide_ind], faceI)
    {
        scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
        scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
        /// 1
        //gTrue[faceI] = - 100000 * faceX * faceZ;
        /// 2
        //if (faceX < moldWidth / 2)
        //{
        //    gTrue[faceI] = - 100000 * faceX * faceZ;
        //}
        //else
        //{
        //    gTrue[faceI] = - 100000 * (moldWidth / 2) * faceZ;
        //    gTrue[faceI] -= - 100000 * (faceX - (moldWidth / 2)) * faceZ;
        //}
        //gTrue[faceI] -= 10000;
        /// 3
        //gTrue[faceI] = 2680 - 335 * Foam::sqrt((moldHeight - faceZ) / castingSpeed);
        //gTrue[faceI] *= - 1000;
        /// 4
        gTrue[faceI] = 100000 * ((faceX - moldWidth / 2) * (faceX - moldWidth / 2) -
                                 moldWidth * moldWidth / 4) - 50000 * faceZ;
        gTrue[faceI] -= 10000;
    }
    set_valueFraction();
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
    readThermocouples();
    Tmeas = fieldValueAtThermocouples(T);

    // Introducing error in the measurements
    //Tmeas += ITHACAutilities::rand(Tmeas.size(), 1, -2, 2);
    //Eigen::VectorXd measurementsError(Tmeas.size());
    //for(int i = 0; i < Tmeas.size(); i++)
    //{
    //    measurementsError(i) = Tmeas.mean() * 0.02 * stochastic::set_normal_random(0.0, 1.0);
    //}
    //Tmeas += measurementsError;

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


void inverseLaplacianProblem::assignDirectBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(T, patchI, - g / k);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
}

void inverseLaplacianProblem::assignAdjointBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& lambda = _lambda();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(lambda, patchI, homogeneousBCcoldSide, refGrad,
                                           valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(lambda, patchI, homogeneousBC);
        }
        else
        {
            ITHACAutilities::assignBC(lambda, patchI, homogeneousBC);
        }
    }
    ITHACAutilities::assignIF(lambda, homogeneousBC);
}

volScalarField inverseLaplacianProblem::assignAdjointBCandSource()
{
    volScalarField& lambda = _lambda();
    assignAdjointBC();
    dimensionedScalar sourceDim("sourceDim", dimensionSet(1, -1, -3, -1, 0, 0, 0),
                                1);
    autoPtr<volScalarField> f_
    (
        new volScalarField("f", lambda)
    );
    volScalarField& f = f_();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            f.ref()[interpolationPlane.cellID [cellI]] =
                interpolationPlane.Tdiff[cellI] * k;
        }
    }
    else
    {
        for (int i = 0; i < thermocouplesCellID.size(); i++)
        {
            f.ref()[thermocouplesCellID [i]] = Tdiff(i) * k;
        }
    }

    return f * sourceDim;
}

void inverseLaplacianProblem::assignSensitivityBC()
{
    fvMesh& mesh = _mesh();
    volScalarField& deltaT = _deltaT();
    set_valueFraction();
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(deltaT, patchI, homogeneousBCcoldSide, refGrad,
                                           valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(deltaT, patchI, - P / k);
        }
        else
        {
            ITHACAutilities::assignBC(deltaT, patchI, homogeneousBC);
        }
    }
}

void inverseLaplacianProblem::solveAdditional()
{
    restart();
    Tad_base.resize(0);
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    volScalarField Tad(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = - Tf;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
        }
        else
        {
            ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
        }
    }
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
                fvm::laplacian(DT, Tad)
            );
            TEqn.solve();
        }
    }

    addSol = fieldValueAtThermocouples(Tad);
    Tad_base.append(Tad);
    ITHACAstream::exportSolution(Tad, "1",
                                 folderOffline,
                                 "Tad");
}

void inverseLaplacianProblem::solveDirect()
{
    restart();
    assignDirectBC();
    solve("direct");
}


void inverseLaplacianProblem::solveAdjoint()
{
    restart();
    volScalarField& lambda = _lambda();
    Foam::Time& runTime = _runTime();
    volScalarField f = assignAdjointBCandSource();
    simpleControl& simple = _simple();
#if OFVER == 6

    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        //Info << "Time = " << runTime.timeName() << nl << endl;
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::laplacian(DT, lambda) ==  - f
            );
            TEqn.solve();
        }
    }
}


void inverseLaplacianProblem::solveSensitivity()
{
    restart();
    assignSensitivityBC();
    solve("sensitivity");
}

void inverseLaplacianProblem::solve(const char* problemID)
{
    volScalarField& T = _T();
    volScalarField& deltaT = _deltaT();
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();

    if (strcmp( problemID, "direct") == 0)
    {
        ITHACAutilities::assignIF(T, homogeneousBC);
    }
    else if (strcmp( problemID, "sensitivity") == 0)
    {
        ITHACAutilities::assignIF(deltaT, homogeneousBC);
    }
    else
    {
        Info << "Problem name should be direct or sensitivity" << endl;
        exit(10);
    }

#if OFVER == 6

    while (simple.loop(runTime))
#else
    while (simple.loop())
#endif
    {
        while (simple.correctNonOrthogonal())
        {
            if (strcmp( problemID, "direct") == 0)
            {
                fvScalarMatrix TEqn
                (
                    fvm::laplacian(DT, T)
                );
                TEqn.solve();
            }
            else if (strcmp( problemID, "sensitivity") == 0)
            {
                fvScalarMatrix TEqn
                (
                    fvm::laplacian(DT, deltaT)
                );
                TEqn.solve();
            }
        }
    }
}

void inverseLaplacianProblem::readThermocouples()
{
    if (!thermocouplesRead)
    {
        word fileName = "./thermocouplesCellsID";

        if (ITHACAutilities::check_file(fileName + "_mat.txt"))
        {
            Info << "Reading thermocouples cells from file" << endl;
            Eigen::MatrixXi TCmatrix = ITHACAstream::readMatrix(fileName + "_mat.txt").cast
                                       <int> ();
            thermocouplesCellID = Foam2Eigen::EigenMatrix2List(TCmatrix);
        }
        else
        {
            Info << "Defining positions of thermocouples" << endl;
            fvMesh& mesh = _mesh();
            volScalarField& T = _T();
            thermocouplesCellID.resize(thermocouplesPos.size());
            forAll(thermocouplesPos, tcI)
            {
                thermocouplesCellID[tcI] = mesh.findCell(thermocouplesPos[tcI]);
            }
            volScalarField thermocouplesField(T);
            ITHACAutilities::assignIF(thermocouplesField, homogeneousBC);
            forAll(thermocouplesCellID, tcI)
            {
                thermocouplesField.ref()[thermocouplesCellID[tcI]] = 1;
            }
            ITHACAstream::exportSolution(thermocouplesField, "1", "./ITHACAoutput/debug/",
                                         "thermocouplesField,");
            Eigen::MatrixXi thermocouplesCellID_eigen = Foam2Eigen::List2EigenMatrix(
                        thermocouplesCellID);
            ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
                                       "eigen", "./");
        }

        thermocouplesRead = 1;
	thermocouplesNum = thermocouplesPos.size();
    }
    else
    {
        WarningInFunction << "readThermocouples function called twice." << endl;
        WarningInFunction << "I am not doing the second reading." << endl;
    }
}

Eigen::VectorXd inverseLaplacianProblem::fieldValueAtThermocouples(
    volScalarField& field)
{
    if (!thermocouplesRead)
    {
        readThermocouples();
    }

    fvMesh& mesh = _mesh();
    dictionary interpolationDict =
        mesh.solutionDict().subDict("interpolationSchemes");
    autoPtr<Foam::interpolation<scalar>> fieldInterp =
                                          Foam::interpolation<scalar>::New(interpolationDict, field);
    Eigen::VectorXd fieldInt;
    fieldInt.resize(thermocouplesPos.size());
    forAll(thermocouplesPos, tcI)
    {
        //label cellI = mesh.findCell(thermocouplesPos[tcI]);
        fieldInt(tcI) = fieldInterp->interpolate(thermocouplesPos[tcI],
                        thermocouplesCellID[tcI]);
    }
    return fieldInt;
}


//void inverseLaplacianProblem::readThermocouples()
//{
//    fvMesh& mesh = _mesh();
//
//    if (!thermocouplesRead)
//    {
//
//  //Define the locations of the thermocouples
//        //All dimensions are in [m]
//
//  //Analytical benchmark
//  label thermocouplesRows = 4;
//  label thermocouplesCols = 4;
//  Eigen::VectorXd thermocouplesX = Eigen::ArrayXd::LinSpaced(thermocouplesCols, 0.2, 0.8);
//  Eigen::VectorXd thermocouplesZ = Eigen::ArrayXd::LinSpaced(thermocouplesRows, 0.2, 0.8);
//  scalar thermocouplesY = 0.2;
//  thermocouplesNum = thermocouplesCols * thermocouplesRows;
//
//  /*
//  // Danieli's thermocouples
//        label thermocouplesRows = 10;
//        label thermocouplesCols = 9;
//        Eigen::VectorXd thermocouplesX = Eigen::ArrayXd::LinSpaced(9, 0.197, 1.773);
//        Eigen::VectorXd thermocouplesZ;
//        thermocouplesZ.resize(thermocouplesRows);
//        thermocouplesZ << 0.120, 0.2215, 0.3225, 0.4335, 0.5445, 0.6555, 0.7665, 0.8775,
//                       0.9885, 1.08;
//        scalar thermocouplesY = 0.015;
//        */
//
//  thermocouplesCellID.resize(thermocouplesNum);
//        thermocouplesCellProc.resize(thermocouplesNum);
//        thermocouplesCellC.resize(thermocouplesNum);
//        Tmeas.resize(thermocouplesNum);
//        Tdirect.resize(thermocouplesNum);
//        Tsens.resize(thermocouplesNum);
//        Tdiff.resize(thermocouplesNum);
//        label cellID = 0;
//
//        for (label rowI = 0; rowI < thermocouplesRows; rowI++)
//        {
//            for (label colI = 0; colI < thermocouplesCols; colI++)
//            {
//                thermocouplesCellID[cellID] =
//                    mesh.findCell(point(thermocouplesX(colI),
//                                        thermocouplesY,
//                                        thermocouplesZ(rowI)));
//                cellID++;
//            }
//        }
//
//        //thermocouplesCellID [0]  = mesh.findCell(point( 0.3, 0.05, 0.9 ));
//        //thermocouplesCellID [1]  = mesh.findCell(point( 0.6, 0.05, 0.9 ));
//        //thermocouplesCellID [2]  = mesh.findCell(point( 0.9, 0.05, 0.9 ));
//        //thermocouplesCellID [3]  = mesh.findCell(point( 1.2, 0.05, 0.9 ));
//        //thermocouplesCellID [4]  = mesh.findCell(point( 1.5, 0.05, 0.9 ));
//        //thermocouplesCellID [5]  = mesh.findCell(point( 0.3, 0.05, 0.6 ));
//        //thermocouplesCellID [6]  = mesh.findCell(point( 0.6, 0.05, 0.6 ));
//        //thermocouplesCellID [7]  = mesh.findCell(point( 0.9, 0.05, 0.6 ));
//        //thermocouplesCellID [8]  = mesh.findCell(point( 1.2, 0.05, 0.6 ));
//        //thermocouplesCellID [9]  = mesh.findCell(point( 1.5, 0.05, 0.6 ));
//        //thermocouplesCellID [10] = mesh.findCell(point( 0.3, 0.05, 0.3 ));
//        //thermocouplesCellID [11] = mesh.findCell(point( 0.6, 0.05, 0.3 ));
//        //thermocouplesCellID [12] = mesh.findCell(point( 0.9, 0.05, 0.3 ));
//        //thermocouplesCellID [13] = mesh.findCell(point( 1.2, 0.05, 0.3 ));
//        //thermocouplesCellID [14] = mesh.findCell(point( 1.5, 0.05, 0.3 ));
//        //Define thermocouples plane
//        unsigned firstCell = 1;
//        forAll(thermocouplesCellID, cellI)
//        {
//            if (thermocouplesCellID[cellI] != -1)
//            {
//                thermocouplesCellProc[cellI] = Pstream::myProcNo();
//
//                if (firstCell)
//                {
//                    interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
//                    interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0);
//                    interpolationPlane.Y = mesh.C()[thermocouplesCellID[cellI]].component(1);
//                    interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
//                    interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
//                    firstCell = 0;
//                }
//
//                if (mesh.C()[thermocouplesCellID[cellI]].component(0) < interpolationPlane.minX)
//                {
//                    interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
//                }
//
//                if (mesh.C()[thermocouplesCellID[cellI]].component(0) > interpolationPlane.maxX)
//                {
//                    interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0) ;
//                }
//
//                if (mesh.C()[thermocouplesCellID[cellI]].component(2) < interpolationPlane.minZ)
//                {
//                    interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
//                }
//
//                if (mesh.C()[thermocouplesCellID[cellI]].component(2) > interpolationPlane.maxZ)
//                {
//                    interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2) ;
//                }
//            }
//            else
//            {
//                Tmeas (cellI) = 0;
//                thermocouplesCellProc[cellI] = -1;
//            }
//        }
//        /*
//        reduce(interpolationPlane.minX,minOp<double>());
//        reduce(interpolationPlane.minZ,minOp<double>());
//        reduce(interpolationPlane.maxX,maxOp<double>());
//        reduce(interpolationPlane.maxZ,maxOp<double>());
//        if(Pstream:: parRun()) //Look up if multiple processors found the same thermocouple
//        {
//            label proc = 0;
//            label cell = 0;
//            if (Pstream::myProcNo() == 0)
//                   {
//                       List<List<int>> temp(nProcs, List<int>(thermocouplesNum));
//                       temp[0] = thermocouplesCellProc;
//                       for(label i=1; i<nProcs; i++)
//                       {
//                           // create the input stream from processor i
//                           IPstream vStream(Pstream::commsTypes::blocking, i);
//                           vStream >> temp[i];
//                       }
//                label flag;
//                for(label i=0; i<thermocouplesNum; i++)
//                {
//                    flag = 0;
//                    //TODO: proc and cell have to become vectors in case of multiple finds
//               for(label j=0; j<nProcs; j++)
//                   {
//                      if(temp[j][i] != -1 && flag == 0)
//                      {
//                          flag = 1;
//                      }
//                      else if(temp[j][i] != -1 && flag == 1)
//                      {
//                          proc = j;
//                          cell = i;
//                          Info << "Double thermocouple cell found." <<
//                       "Removing Proc = " << proc<< " Thermocouple = "<<  cell<< endl;
//                      }
//                   }
//                }
//
//                   }
//                   else
//                   {
//                       // create the stream to send to the main proc
//                       OPstream vectorStream
//                       (
//                           Pstream::commsTypes::blocking, 0
//                       );
//                       vectorStream << thermocouplesCellProc;
//                   }
//            reduce(proc, sumOp<label>());
//            reduce(cell, sumOp<label>());
//
//            if(Pstream::myProcNo() == proc)
//            {
//                thermocouplesCellID[cell] = -1;
//                Tmeas (cell) = 0;
//            }
//            if(interpolation)
//            {
//                forAll(thermocouplesCellID, cellI)
//                {
//                    reduce(Tmeas (cellI), sumOp<double>());
//                       }
//            }
//        }
//        */
//    }
//    else
//    {
//        WarningInFunction << "readThermocouples function called twice." << endl;
//        WarningInFunction << "I am not doing the second reading." << endl;
//    }
//}

void inverseLaplacianProblem::defineThermocouplesPlane()
{
    Info << "Defining the plane for measurements interpolation" << endl;
    fvMesh& mesh = _mesh();
    //Define thermocouples plane
    bool firstCell = 1;
    forAll(thermocouplesCellID, cellI)
    {
        if (thermocouplesCellID[cellI] != -1)
        {
            thermocouplesCellProc[cellI] = Pstream::myProcNo();

            if (firstCell)
            {
                interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
                interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0);
                interpolationPlane.Y = mesh.C()[thermocouplesCellID[cellI]].component(1);
                interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
                interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
                firstCell = 0;
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(0) < interpolationPlane.minX)
            {
                interpolationPlane.minX = mesh.C()[thermocouplesCellID[cellI]].component(0);
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(0) > interpolationPlane.maxX)
            {
                interpolationPlane.maxX = mesh.C()[thermocouplesCellID[cellI]].component(0) ;
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(2) < interpolationPlane.minZ)
            {
                interpolationPlane.minZ = mesh.C()[thermocouplesCellID[cellI]].component(2);
            }

            if (mesh.C()[thermocouplesCellID[cellI]].component(2) > interpolationPlane.maxZ)
            {
                interpolationPlane.maxZ = mesh.C()[thermocouplesCellID[cellI]].component(2) ;
            }
        }
        else
        {
            Tmeas (cellI) = 0;
            thermocouplesCellProc[cellI] = -1;
        }
    }
}

void inverseLaplacianProblem::differenceBetweenDirectAndMeasure()
{
    volScalarField& T = _T();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            interpolationPlane.Tdirect[cellI] =
                T.internalField()[interpolationPlane.cellID [cellI]];
        }
        interpolationPlane.Tdiff = interpolationPlane.Tdirect -
                                   interpolationPlane.Tmeas;
    }
    else
    {
        Tdirect = fieldValueAtThermocouples(T);
        Tdiff = Tdirect - Tmeas;
    }
}

void inverseLaplacianProblem::sensibilitySolAtThermocouplesLocations()
{
    volScalarField& deltaT = _deltaT();

    if (interpolation)
    {
        forAll(interpolationPlane.cellID, cellI)
        {
            interpolationPlane.Tsens[cellI] =
                deltaT.internalField()[interpolationPlane.cellID [cellI]];
        }
    }
    else
    {
        Tsens = fieldValueAtThermocouples(deltaT);
    }
}

int inverseLaplacianProblem::conjugateGradient()
{
    set_g();
    set_valueFraction();
    cgIter = 0;
    J = 0;
    P = g;
    gradJ = g;       //Gradient of the cost function [W/m2]
    gamma = 0.0;
    gamma_den = 0.0;
    label sampleI = 1;
    gList.resize(0);

    while (cgIter < cgIterMax)
    {
        Info << "Iteration " << cgIter + 1 << endl;
        restart();
        solveDirect();

        if (saveSolInLists && cgIter == 0)
        {
            gList.append(g);
        }

        volScalarField& T = _T();
        ITHACAstream::exportSolution(T, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", T.name());
        differenceBetweenDirectAndMeasure();

        if (conjugateGradientConvergenceCheck())
        {
            Jlist.conservativeResize(cgIter + 1, 1);
            Jlist(cgIter) = J;
            ITHACAstream::exportMatrix(Jlist, "costFunctionFull", "eigen", "./");
            return (1);
        }

        Jlist.conservativeResize(cgIter + 1, 1);
        Jlist(cgIter) = J;
        solveAdjoint();
        volScalarField& lambda = _lambda();
        ITHACAstream::exportSolution(lambda, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", lambda.name());
        computeGradJ();
        searchDirection();
        solveSensitivity();
        volScalarField& deltaT = _deltaT();
        ITHACAstream::exportSolution(deltaT, std::to_string(sampleI),
                                     "./ITHACAoutput/CGtest/", deltaT.name());
        sensibilitySolAtThermocouplesLocations();
        computeSearchStep();
        updateHeatFlux();
        // I save the fields at each iteration in the Offline folder
        // the counter should not overload the solution
        offlineSolutionI++;

        // I the online phase I only save a solution each 5 iterations and the last one
        if (offlinePhase && offlineSolutionI % 5 == 0 )
        {
            label saveIndex = offlineSolutionI / 5;
            writeFields(saveIndex, "./ITHACAoutput/Offline/");
        }

        if (saveSolInLists)
        {
            volScalarField& T = _T();
            volScalarField& lambda = _lambda();
            volScalarField& deltaT = _deltaT();
            gList.append(g);
            Tfield.append(T);
            lambdaField.append(lambda);
            deltaTfield.append(deltaT);
            sampleI++;
        }

        cgIter++;
    }

    return (0);
}

void inverseLaplacianProblem::computeGradJ()
{
    fvMesh& mesh = _mesh();
    volScalarField& lambda = _lambda();
    gradJ_L2norm = 0;
    forAll (lambda.boundaryField()[hotSide_ind], faceI)
    {
        gradJ [faceI] = - lambda.boundaryField()[hotSide_ind][faceI];
        gradJ_L2norm += gradJ[faceI] * gradJ[faceI]  *
                        mesh.magSf().boundaryField()[hotSide_ind][faceI];
    }
    gradJ_L2norm = Foam::sqrt(gradJ_L2norm);
    Info << "gradJ L2norm = " << gradJ_L2norm << endl;
}

void inverseLaplacianProblem::searchDirection()
{
    fvMesh& mesh = _mesh();
    gamma = 0.0;
    scalar gammaNum = 0;
    forAll (mesh.magSf().boundaryField()[hotSide_ind], faceI)
    {
        gammaNum += gradJ [faceI] * gradJ [faceI] *
                    mesh.magSf().boundaryField()[hotSide_ind][faceI];
    }

    if (cgIter > 0)
    {
        //reduce(gamma, sumOp<double>());
        gamma = gammaNum / gamma_den;
        Info << "gamma = " << gamma << endl;
    }

    P = gradJ + gamma * P; //Updating P
    gamma_den = gammaNum;
    //reduce(gamma_den, sumOp<double>());
}

void inverseLaplacianProblem::computeSearchStep()
{
    if (interpolation)
    {
        List<scalar> temp = interpolationPlane.Tdiff *
                            interpolationPlane.Tsens;
        beta = 0.0;
        forAll(interpolationPlane.cellVol, cellI)
        {
            beta += interpolationPlane.cellVol [cellI] * temp [cellI];
        }
        //reduce(beta, sumOp<double>());
        temp = interpolationPlane.Tsens * interpolationPlane.Tsens;
        scalar betaDiv = 0.0;
        forAll(interpolationPlane.cellVol, cellI)
        {
            betaDiv += interpolationPlane.cellVol [cellI] * temp [cellI];
        }
        //reduce(betaDiv, sumOp<double>());
        beta = beta / betaDiv;
        temp.clear();
    }
    else
    {
        beta = Tdiff.dot(Tsens);
        //reduce(beta, sumOp<double>());
        double betaDiv = Tsens.dot(Tsens);
        //reduce(betaDiv, sumOp<double>());
        beta = beta / betaDiv;
    }

    Info << "beta = " << beta << endl;
}

void inverseLaplacianProblem::updateHeatFlux()
{
    g = g - beta * P;
}


int inverseLaplacianProblem::conjugateGradientConvergenceCheck()
{
    double Jold = J;

    if (interpolation)
    {
        List<scalar> sqTdiff;
        sqTdiff = interpolationPlane.Tdiff * interpolationPlane.Tdiff;
        J = 0.0;
        forAll(sqTdiff, cellI)
        {
            J += 0.5 * sqTdiff[cellI] * interpolationPlane.cellVol[cellI];
        }
        sqTdiff.clear();
    }
    else
    {
        J = 0.5 * Tdiff.dot(Tdiff);
    }

    //reduce(J, sumOp<double>());
    Info << "J = " << J << endl;

    if (J <= Jtol)
    {
        Info << "Convergence reached in " << cgIter << " iterations" << endl;
        return (1);
    }
    //else if (cgIter > 0 && gradJ_L2norm <= Jtol)
    //{
    //    Info << "Stopping criteria on the gradient of the cost functional" <<
    //         endl << "met in " << cgIter << " iterations" << endl;
    //    return (1);
    //}
    else if (Foam::mag((Jold - J) / J) <= JtolRel)
    {
        Info << "Relative tolerance criteria meet in " << cgIter << " iterations" <<
             endl;
        Info << "|Jold - J| / |J| = " << Foam::mag((Jold - J) / J) << endl;
        return (1);
    }
    else
    {
        return (0);
    }
}



int inverseLaplacianProblem::isInPlane(double cx, double cy, double cz,
                                       Foam::vector thermocoupleCellDim)
{
    return (cx >= interpolationPlane.minX - thermocoupleCellDim[0] / 4 &&
            cy >= interpolationPlane.Y - thermocoupleCellDim[1] / 4 &&
            cy <= interpolationPlane.Y + thermocoupleCellDim[1] / 4 &&
            cz >= interpolationPlane.minZ - thermocoupleCellDim[2] / 4 &&
            cx <= interpolationPlane.maxX + thermocoupleCellDim[0] / 4 &&
            cz <= interpolationPlane.maxZ + thermocoupleCellDim[2] / 4
           );
}

void inverseLaplacianProblem::writeFields(label folderNumber,
        const char* folder)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    volScalarField& lambda = _lambda();
    volScalarField& deltaT = _deltaT();
    autoPtr<volScalarField> gVolField_
    (
        new volScalarField("g", T)
    );
    volScalarField& gVolField = gVolField_();
    ITHACAutilities::assignIF(gVolField, homogeneousBC);
    //Access the mesh information for the boundary
    const polyPatch& cPatch = mesh.boundaryMesh()[hotSide_ind];
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    forAll(cPatch, faceI)
    {
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        gVolField[faceOwner] = g[faceI];
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


Foam::vector inverseLaplacianProblem::cellDim(const faceList& ff,
        const pointField& pp,
        const cell& cc, labelList pLabels, pointField pLocal)
{
    forAll (pLabels, pointi)
    pLocal[pointi] = pp[pLabels[pointi]];
    double  xDim = Foam::max(pLocal & Foam::vector(1, 0, 0))
                   - Foam::min(pLocal & Foam::vector(1, 0, 0));
    double  yDim = Foam::max(pLocal & Foam::vector(0, 1, 0))
                   - Foam::min(pLocal & Foam::vector(0, 1, 0));
    double  zDim = Foam::max(pLocal & Foam::vector(0, 0, 1))
                   - Foam::min(pLocal & Foam::vector(0, 0, 1));
    Foam::vector dim (xDim, yDim, zDim);
    return dim;
}

//Interpolates the values in Tmeas on the interpolation plane defined in readThermocouples()
void inverseLaplacianProblem::thermocouplesInterpolation()
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    DataTable thermocouplesSamples;
    DenseVector x(2);
    unsigned i = 0;

    if (!interpolationPlaneDefined)
    {
        interpolationPlane.thermocoupleX.resize(thermocouplesNum);
        interpolationPlane.thermocoupleZ.resize(thermocouplesNum);

        //find the dimensions of the cells containing the thermocouples
        //I am assuming all thermocouples are a the same y coordinate
        //I am assuming all cells have the same dimensions
        if ( Pstream::master() == true )
        {
            unsigned flag = 1;

            while (flag == 1)
            {
                if (thermocouplesCellID [i] != -1)
                {
                    labelList pLabels(mesh.cells()[thermocouplesCellID [i]].labels(mesh.faces()));
                    pointField pLocal(pLabels.size(), Foam::vector::zero);
                    interpolationPlane.thermocoupleCellDim = cellDim (mesh.faces(),
                            mesh.points(),
                            mesh.cells()[thermocouplesCellID [i]],
                            pLabels,
                            pLocal);
                    flag = 0;
                }

                i++;
            }
        }

        //reduce(interpolationPlane.thermocoupleCellDim, sumOp<vector>());
        forAll(thermocouplesCellID, thermocoupleI)
        {
            if (thermocouplesCellID[thermocoupleI] != -1)
            {
                interpolationPlane.thermocoupleX [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(0);
                interpolationPlane.thermocoupleZ [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(2);
            }
            else
            {
                interpolationPlane.thermocoupleX [thermocoupleI] = 0;
                interpolationPlane.thermocoupleZ [thermocoupleI] = 0;
            }
        }
        //reduce(interpolationPlane.thermocoupleX, sumOp<List<scalar>>());
        //reduce(interpolationPlane.thermocoupleZ, sumOp<List<scalar>>());
    }

    forAll(thermocouplesCellID, thermocoupleI)
    {
        x(0) = interpolationPlane.thermocoupleX[thermocoupleI];
        x(1) = interpolationPlane.thermocoupleZ[thermocoupleI];
        thermocouplesSamples.addSample(x, Tmeas(thermocoupleI));
    }
    std::cout << Tmeas << std::endl;
    RBFSpline rbfspline(thermocouplesSamples, RadialBasisFunctionType::GAUSSIAN);
    auto inPlaneCellID = 0;
    forAll(T.internalField(), cellI)
    {
        auto cx = mesh.C()[cellI].component(Foam::vector::X);
        auto cy = mesh.C()[cellI].component(Foam::vector::Y);
        auto cz = mesh.C()[cellI].component(Foam::vector::Z);

        if (!interpolationPlaneDefined)
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                auto planeSize = interpolationPlane.cellID.size() + 1;
                interpolationPlane.cellID.resize (planeSize);
                interpolationPlane.Tmeas.resize  (planeSize);
                interpolationPlane.Tdirect.resize(planeSize);
                interpolationPlane.Tdiff.resize  (planeSize);
                interpolationPlane.Tsens.resize  (planeSize);
                interpolationPlane.cellVol.resize(planeSize);
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.cellID [planeSize - 1] = cellI;
                interpolationPlane.Tmeas  [planeSize - 1] =
                    rbfspline.eval(x);
                interpolationPlane.cellVol[planeSize - 1] =
                    mesh.V()[cellI];
            }
        }
        else
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.Tmeas  [inPlaneCellID] =
                    rbfspline.eval(x);
                inPlaneCellID++;
            }
        }
    }
    interpolationPlaneDefined = 1;
}

//Interpolates the values in Tmeas on the interpolation plane defined in readThermocouples()
void inverseLaplacianProblem::thermocouplesInterpolation(
    DenseMatrix& RBFweights, DenseMatrix& RBFbasis)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    DataTable thermocouplesSamples;
    DenseVector x(2);
    unsigned i = 0;

    if (!interpolationPlaneDefined)
    {
        interpolationPlane.thermocoupleX.resize(thermocouplesNum);
        interpolationPlane.thermocoupleZ.resize(thermocouplesNum);

        //find the dimensions of the cells containing the thermocouples
        //I am assuming all thermocouples are a the same y coordinate
        //I am assuming all cells have the same dimensions
        if ( Pstream::master() == true )
        {
            unsigned flag = 1;

            while (flag == 1)
            {
                if (thermocouplesCellID [i] != -1)
                {
                    labelList pLabels(mesh.cells()[thermocouplesCellID [i]].labels(mesh.faces()));
                    pointField pLocal(pLabels.size(), Foam::vector::zero);
                    interpolationPlane.thermocoupleCellDim = cellDim (mesh.faces(),
                            mesh.points(),
                            mesh.cells()[thermocouplesCellID [i]],
                            pLabels,
                            pLocal);
                    flag = 0;
                }

                i++;
            }
        }

        //reduce(interpolationPlane.thermocoupleCellDim, sumOp<vector>());
        forAll(thermocouplesCellID, thermocoupleI)
        {
            if (thermocouplesCellID[thermocoupleI] != -1)
            {
                interpolationPlane.thermocoupleX [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(0);
                interpolationPlane.thermocoupleZ [thermocoupleI] =
                    mesh.C()[thermocouplesCellID [thermocoupleI]].component(2);
            }
            else
            {
                interpolationPlane.thermocoupleX [thermocoupleI] = 0;
                interpolationPlane.thermocoupleZ [thermocoupleI] = 0;
            }
        }
        //reduce(interpolationPlane.thermocoupleX, sumOp<List<scalar>>());
        //reduce(interpolationPlane.thermocoupleZ, sumOp<List<scalar>>());
    }

    forAll(thermocouplesCellID, thermocoupleI)
    {
        x(0) = interpolationPlane.thermocoupleX[thermocoupleI];
        x(1) = interpolationPlane.thermocoupleZ[thermocoupleI];
        thermocouplesSamples.addSample(x, Tmeas(thermocoupleI));
    }
    std::cout << Tmeas << std::endl;
    RBFSpline rbfspline(thermocouplesSamples, RadialBasisFunctionType::GAUSSIAN);
    auto inPlaneCellID = 0;
    forAll(T.internalField(), cellI)
    {
        auto cx = mesh.C()[cellI].component(Foam::vector::X);
        auto cy = mesh.C()[cellI].component(Foam::vector::Y);
        auto cz = mesh.C()[cellI].component(Foam::vector::Z);

        if (!interpolationPlaneDefined)
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                auto planeSize = interpolationPlane.cellID.size() + 1;
                interpolationPlane.cellID.resize (planeSize);
                interpolationPlane.Tmeas.resize  (planeSize);
                interpolationPlane.Tdirect.resize(planeSize);
                interpolationPlane.Tdiff.resize  (planeSize);
                interpolationPlane.Tsens.resize  (planeSize);
                interpolationPlane.cellVol.resize(planeSize);
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.cellID [planeSize - 1] = cellI;
                interpolationPlane.Tmeas  [planeSize - 1] =
                    rbfspline.eval(x);
                interpolationPlane.cellVol[planeSize - 1] =
                    mesh.V()[cellI];
            }
        }
        else
        {
            if (isInPlane(cx, cy, cz, interpolationPlane.thermocoupleCellDim))
            {
                x(0) = cx;
                x(1) = cz;
                interpolationPlane.Tmeas  [inPlaneCellID] =
                    rbfspline.eval(x);
                inPlaneCellID++;
            }
        }
    }
    interpolationPlaneDefined = 1;
}

void inverseLaplacianProblem::offlineSolve()
{
    volScalarField& T = _T();
    volScalarField& lambda = _lambda();
    volScalarField& deltaT = _deltaT();

    if (offline)
    {
        ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline/");
        ITHACAstream::read_fields(lambdaField, lambda, "./ITHACAoutput/Offline/");
        ITHACAstream::read_fields(deltaTfield, deltaT, "./ITHACAoutput/Offline/");
    }
    else
    {
        offlinePhase = 1;
        set_g();
        readThermocouples();
        set_valueFraction();
        saveSolInLists = 0;
        offlineSolutionI = 1;

        for (label i = 0; i < muSamples.cols(); i++)
        {
            Info << "Sample number " << i + 1 << endl;
            Tmeas = muSamples.col(i);

            if (interpolation)
            {
                Info << "Interpolating thermocouples measurements in the " <<
                     "plane defined by them" << endl;
                thermocouplesInterpolation();
            }
            else
            {
                Info << "NOT interpolating thermocouples measurements" << endl;
            }

            if (!conjugateGradient())
            {
                Info << "Conjugate gradient method did not converge" <<
                     "Exiting" << endl;
                exit(10);
            }
        }

        saveSolInLists = 0;
        offlinePhase = 0;
    }
}

void inverseLaplacianProblem::samplingHeatFluxMarquardt(word folder)
{
    Eigen::MatrixXd gWeightsOffline;

    if (ITHACAutilities::check_folder(folder))
    {
        Info << "Levemberg-Marquardt offline already computed" << endl;
        MarquardtHeatFluxSampling = 1;
        //ITHACAstream::read_fields(Tfield, T, folder);
        set_gParametrized("rbf");
        Eigen::MatrixXd gBaseFunctionsEigen;
        gBaseFunctionsEigen.resize(gBaseFunctions[0].size(), gBaseFunctions.size());
        Info << "gBaseFunctionsEigen size = " << gBaseFunctionsEigen.rows() << "," <<
             gBaseFunctionsEigen.cols() << endl;

        for (label i = 0; i < gBaseFunctions.size() ; i++)
        {
            for (label j = 0; j < gBaseFunctions[i].size() ; j++)
            {
                gBaseFunctionsEigen(j, i) = gBaseFunctions[i][j];
            }
        }

        gWeightsOffline = ITHACAstream::readMatrix(folder + "/gWeightsOffline_mat.txt");
    }
    else
    {
        Info << "Sampling Tmeas Levemberg-Marquardt" << endl;
        offlinePhase = 1;
        //set_g();
        readThermocouples();
        set_valueFraction();
        saveSolInLists = 0;
        offlineSolutionI = 1;
        label updateJacobian = 0;
        set_gParametrized("rbf");
        MarquardtMethodSetUp();

        for (label i = 0; i < muSamples.cols(); i++)
        {
            Info << "Marquardt offline, Sample number " << i + 1 << endl;
            Tmeas = muSamples.col(i);
            std::cout << "Tmeas = " << Tmeas << endl;
            MarquardtMethod(1, folder, updateJacobian);
            gWeightsOffline.conservativeResize(gWeights.size(), i + 1);

            for (label j = 0; j < gWeights.size(); j++)
            {
                gWeightsOffline(j, i) = gWeights[j];
            }
        }

        ITHACAstream::exportMatrix(gWeightsOffline, "gWeightsOffline", "eigen", folder);
        saveSolInLists = 0;
        offlinePhase = 0;
        MarquardtHeatFluxSampling = 1;
    }
}

void inverseLaplacianProblem::heatFluxPodMarquardt(word folder)
{
    Eigen::MatrixXd gWeightsSamples;
    M_Assert(MarquardtHeatFluxSampling,
             "Before calling podMarquardt you must call offlineSolveMarquardt.");

    if (MarquardtHeatFluxPOD) //ITHACAutilities::check_folder(folder))
    {
        Info << "Levemberg-Marquardt POD of the heat flux already computed" << endl;
        // COMPLETE
    }
    else
    {
        set_gParametrized("rbf");
        Eigen::MatrixXd gBaseFunctionsEigen;
        gBaseFunctionsEigen.resize(gBaseFunctions[0].size(), gBaseFunctions.size());
        Info << "gBaseFunctionsEigen size = " << gBaseFunctionsEigen.rows() << "," <<
             gBaseFunctionsEigen.cols() << endl;

        for (label i = 0; i < gBaseFunctions.size() ; i++)
        {
            for (label j = 0; j < gBaseFunctions[i].size() ; j++)
            {
                gBaseFunctionsEigen(j, i) = gBaseFunctions[i][j];
            }
        }

        gWeightsSamples =
            ITHACAstream::readMatrix("./ITHACAoutput/offlineMarquardt/gWeightsOffline_mat.txt");
        Info << "gWeightsSamples size = " << gWeightsSamples.rows() << "," <<
             gWeightsSamples.cols() << endl;
        Eigen::MatrixXd gSamples = gBaseFunctionsEigen * gWeightsSamples;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(gSamples, Eigen::ComputeThinU);
        Eigen::MatrixXd singularValues = svd.singularValues();
        Eigen::MatrixXd U = svd.matrixU();
        Info << "U size = " << U.rows() << "," << U.cols() << endl;
        ITHACAstream::exportMatrix(singularValues, "singularValues", "eigen", folder);
        ITHACAstream::exportMatrix(U, "PODbases", "eigen", folder);
        MarquardtHeatFluxPOD = 1;
        label PODout = 100;

        for (label i = 0; i < PODout; i++) //Only save PODout modes
        {
            volScalarField gBase = eigen2Field(U.col(i));
            ITHACAstream::exportSolution(gBase, std::to_string(i + 1), folder,
                                         "gPOD");
        }

        label Nbasis = 10;
        Eigen::MatrixXd gReducedBases = U.leftCols(Nbasis);
        Eigen::MatrixXd gWeightsReduced = gReducedBases.transpose() * gSamples;
        ITHACAstream::exportMatrix(gReducedBases, "gReducedBases", "eigen", folder);
        ITHACAstream::exportMatrix(gWeightsReduced, "gWeightsReduced", "eigen", folder);
        ITHACAstream::exportMatrix(gSamples, "gSamples", "eigen", folder);
    }
}

void inverseLaplacianProblem::MarquardtOffline(word folder)
{
    volScalarField& T = _T();
    Eigen::MatrixXd gWeightsOffline;

    if (ITHACAutilities::check_folder(folder))
    {
        Info << "Levemberg-Marquardt offline already computed" << endl;
    }
    else
    {
        Info << "Computing Levemberg-Marquardt offline" << endl;
        offlinePhase = 1;
        readThermocouples();
        saveSolInLists = 0;
        offlineSolutionI = 1;
        label updateJacobian = 0;
        set_gParametrized("pod");
        MarquardtMethodSetUp();
        Eigen::MatrixXd gWeightsSamples;
        cnpy::load(gWeightsSamples, "./ITHACAoutput/podMarquardt/gWeightSamples.npy");

        for (label i = 0; i < gWeightsSamples.rows(); i++)
        {
            Info << "Sample " << i + 1 << endl;
            forAll(gWeights, weightI)
            {
                gWeights[weightI] = gWeightsSamples(i, weightI);
            }
            update_gParametrized(gWeights);
            Info << "gWeights = " << gWeights << endl;
            solveDirect();
            volScalarField gField = list2Field(g);
            ITHACAstream::exportSolution(gField, std::to_string(i + 1),
                                         folder, "g");
            ITHACAstream::exportSolution(T, std::to_string(i + 1),
                                         folder, "T");
        }
    }

    ITHACAstream::read_fields(Tfield, T, "./ITHACAoutput/Offline_Marquardt/");
}

void inverseLaplacianProblem::restart(word fieldName)
{
    Info << "\nResetting time, mesh and fields: " << fieldName << "\n" << endl;
    _simple.clear();

    if (fieldName == "T" || fieldName == "all")
    {
        _T.clear();
    }

    if (fieldName == "lambda" || fieldName == "all")
    {
        _lambda.clear();
    }

    if (fieldName == "deltaT" || fieldName == "all")
    {
        _deltaT.clear();
    }

    argList& args = _args();
    Time& runTime = _runTime();

    //Reinitializing runTime
    instantList Times = runTime.times();
    runTime.setTime(Times[1], 1);

    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );

    if (fieldName == "T" || fieldName == "all")
    {
        //Info << "ReReading field T\n" << endl;
        _T = autoPtr<volScalarField>
             (
                 new volScalarField
                 (
                     IOobject
                     (
                         "T",
                         runTime.timeName(),
                         mesh,
                         IOobject::MUST_READ,
                         IOobject::AUTO_WRITE
                     ),
                     mesh
                 )
             );
    }

    if (fieldName == "deltaT" || fieldName == "all")
    {
        //Info << "ReReading field deltaT\n" << endl;
        _deltaT = autoPtr<volScalarField>
                  (
                      new volScalarField
                      (
                          IOobject
                          (
                              "deltaT",
                              runTime.timeName(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::AUTO_WRITE
                          ),
                          mesh
                      )
                  );
    }

    if (fieldName == "lambda" || fieldName == "all")
    {
        //Info << "ReReading field lambda\n" << endl;
        _lambda = autoPtr<volScalarField>
                  (
                      new volScalarField
                      (
                          IOobject
                          (
                              "lambda",
                              runTime.timeName(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::AUTO_WRITE
                          ),
                          mesh
                      )
                  );
    }

    Info << "Ready for new computation" << endl;
}

void inverseLaplacianProblem::reconstructT()
{
    Info << "Reconstructing field T" << endl;
    restart();
    update_gParametrized(gWeights);

    volScalarField& T = _T();
    ITHACAutilities::assignIF(T, homogeneousBC);

    forAll(Tbasis, baseI)
    {
        T += gWeights[baseI] * (Tbasis[baseI] + Tad_base[0]);
    }

    T += - Tad_base[0];
}

