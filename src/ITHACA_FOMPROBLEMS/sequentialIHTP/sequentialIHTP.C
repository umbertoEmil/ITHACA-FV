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
/// Source file of the sequentialIHTP class.


#include "sequentialIHTP.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
sequentialIHTP::sequentialIHTP() {}

sequentialIHTP::sequentialIHTP(int argc, char* argv[])
    :
    DT("DT", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1.0)
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
#include "createRegularization.H"
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    startTime = runTime.startTime().value();
    deltaTime = runTime.deltaTValue();
    endTime = runTime.endTime().value();
    Ntimes = (endTime - startTime) / deltaTime;
    timeSteps.resize( Ntimes + 1 );
    forAll(timeSteps, timeI)
    {
        timeSteps[timeI] = startTime + (timeI) * deltaTime;
    }
    //Info << "debug: timeSteps = " << timeSteps << endl;
    //Info << "debug: startTime = " << startTime << endl;
    //Info << "debug: endTime = " << endTime << endl;
    //Info << "debug: deltaTime = " << deltaTime << endl;
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void sequentialIHTP::setDiffusivity(scalar _diff)
{
    diffusivity = _diff;
}

void sequentialIHTP::setSpaceBasis(word type,
        scalar shapeParameter, label Npod)
{
    M_Assert(type == "rbf" || type == "pod", "Only RBF or POD basis implemented");
    if (!thermocouplesRead)
    {
        readThermocouples();
    }
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    Nbasis = thermocouplesNum;
    heatFluxSpaceBasis.resize(Nbasis);

    Info << "\nRadial Basis Functions are used." << endl;
    Info << "The center of each function is at the projection " << endl;
    Info << "of each thermocouple on the boundary hotSide.\n\n";

    int thermocouplesCounter = 0;
    int rbfCenterTimeI = 0;
    scalar maxX =  Foam::max(
            mesh.boundaryMesh()[hotSide_ind].faceCentres().component(Foam::vector::X));
    scalar maxZ =  Foam::max(
            mesh.boundaryMesh()[hotSide_ind].faceCentres().component(Foam::vector::Z));

    forAll(heatFluxSpaceBasis, funcI)
    {
        scalar thermocoupleX =
            mesh.C()[thermocouplesCellID [thermocouplesCounter]].component(0);
        scalar thermocoupleZ =
            mesh.C()[thermocouplesCellID [thermocouplesCounter]].component(2);
        heatFluxSpaceBasis[funcI].resize(T.boundaryField()[hotSide_ind].size());
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            scalar faceX = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
            scalar faceZ = mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();

            scalar radius = Foam::sqrt((faceX - thermocoupleX) * (faceX - 
                    thermocoupleX) / maxX / maxX + (faceZ - thermocoupleZ) * 
                    (faceZ - thermocoupleZ) / maxZ / maxZ);
            heatFluxSpaceBasis[funcI][faceI] = Foam::sqrt(1 + (shapeParameter *
                        radius) * (shapeParameter * radius));

        }
        thermocouplesCounter++;
    }
    if (type == "pod")
    {
        Info << "Using POD space basis" << endl << endl;
        M_Assert(Npod > 0, "Set number of POD basis");
        Nbasis = Npod;
        List<scalar> massVector(T.boundaryField()[hotSide_ind].size());
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            massVector[faceI] = mesh.boundary()[hotSide_ind].magSf()[faceI];
        }
        List<List<scalar>> tempBasis;
        word debugFolder = "./ITHACAoutput/debugParameterizedBasis/";
        ITHACAPOD::getModesSVD(heatFluxSpaceBasis, massVector, tempBasis, Npod, 
                debugFolder);
        forAll(heatFluxSpaceBasis, baseI)
        {
            volScalarField base = list2Field(heatFluxSpaceBasis[baseI], 0.0);
            ITHACAstream::exportSolution(base,
                                         std::to_string(1),
                                         debugFolder,
                                         "RBFbase" + std::to_string(baseI + 1));
        }
        forAll(tempBasis, baseI)
        {
            volScalarField base = list2Field(tempBasis[baseI], 0.0);
            ITHACAstream::exportSolution(base,
                                         std::to_string(1),
                                         debugFolder,
                                         "PODbase" + std::to_string(baseI + 1));
        }
        heatFluxSpaceBasis = tempBasis;
    }
}

void sequentialIHTP::setParametrizedHeatFlux(word spaceBaseFuncType,
        scalar shapeParameter_space, label Npod)
{
    volScalarField& T = _T();
    setSpaceBasis(spaceBaseFuncType, shapeParameter_space, Npod);

    heatFluxTimeBasis.resize(NtimeStepsBetweenSamples);
    forAll(heatFluxTimeBasis, timeI)
    {
        heatFluxTimeBasis[timeI].resize(Nbasis);
    }

    heatFlux.resize(timeSteps.size());
    heatFluxWeights.resize(Nbasis);
    forAll (heatFluxWeights, weigthI)
    {
        heatFluxWeights[weigthI] = 0;
    }
    forAll(timeSteps, timeI)
    {
        heatFlux[timeI].resize(T.boundaryField()[hotSide_ind].size(), 0.0);
    }
}

void sequentialIHTP::updateHeatFlux(List<scalar> weights)
{
    M_Assert(weights.size() == Nbasis,
             "weigths size different from basis functions size");

    label firstTimeI = timeSampleI * NtimeStepsBetweenSamples;
    int lastTimeStep = firstTimeI + NtimeStepsBetweenSamples;

    updateTimeHeatFlux(weights);
    volScalarField& T = _T();
    for(int timeI = firstTimeI; timeI < lastTimeStep; timeI++)
    {    
        forAll (T.boundaryField()[hotSide_ind], faceI)
        {
            heatFlux[timeI + 1][faceI] = 0.0;
            forAll (weights, weightI)
            {
                heatFlux[timeI + 1][faceI] += 
                    heatFluxTimeBasis[timeI - firstTimeI][weightI] * 
                    heatFluxSpaceBasis[weightI][faceI];
            }
        }
        Info << endl << endl;
    }
}

volScalarField sequentialIHTP::list2Field(List<scalar> list,
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

void sequentialIHTP::solveT(volScalarField _initialField, word outputFolder)
{
    Info << "Solving for field T" << endl;
    Ttime.resize(0);
    Info << "\nSolving in the time domain (" << 
        timeSteps[NtimeStepsBetweenSamples * timeSampleI] << ", " << 
        timeSteps[NtimeStepsBetweenSamples + NtimeStepsBetweenSamples * timeSampleI] << 
        "]\n" << endl;
    restartOffline();
    if(timeSampleI == 0)
    {
        ITHACAstream::exportSolution(_initialField, std::to_string(timeSteps[0]),
                                     outputFolder,
                                     "Tsol");
    }
    M_Assert(diffusivity>1e-36, "Set the diffusivity value");
    volScalarField& T = _T();
    ITHACAutilities::assignIF(T, _initialField); 
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();
    fv::options& fvOptions(_fvOptions());
    label timeI = 0;

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        label realTimeStep = timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
        //Info << "debug: realTimeStep = " << realTimeStep << endl;
        assignDirectBC(realTimeStep);

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT * diffusivity, T)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        ITHACAstream::exportSolution(T, std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "Tsol");
        volScalarField gParametrizedField = list2Field(heatFlux[realTimeStep]);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "gSol");
        Ttime.append(T.clone());

        runTime.printExecutionTime(Info);
        runTime.write();
        timeI++;
    }
    Info << "solveT ENDED" << endl;

}

Eigen::VectorXd sequentialIHTP::solveLinSys(List<Eigen::MatrixXd> linSys)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
    //ITHACAregularization::Picard(svd.matrixU(), svd.singularValues(),
    //        linSys[1]);
    M_Assert(linSys.size() == 2, "Linear system has wrong dimension");
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
        Info << "Using TSVD" << endl;
        weigths = ITHACAregularization::TSVD(linSys[0], linSys[1],
                TSVD_filter);
    }
    else if (linSys_solver == "Tikhonov")
    {
        ITHACAregularization::Lcurve(svd.matrixU(), svd.singularValues(),
                linSys[1], linSys_solver, timeSampleI);

        weigths = ITHACAregularization::Tikhonov(svd.matrixU(), svd.singularValues(),
                svd.matrixV(), linSys[1], Tikhonov_filter);
    }
    else if (linSys_solver == "conjugateGradient")
    {
        ITHACAregularization::Lcurve_CG(linSys[0], linSys[1]);
        weigths = ITHACAregularization::conjugateGradient(linSys[0], 
                linSys[1], CG_Nsteps);
    }
    else if (linSys_solver == "PCGLS")
    {
        label Lorder = 1;
        Eigen::MatrixXd W(linSys[0].rows(), Lorder);
        Eigen::MatrixXd L;

        L = ITHACAregularization::get_L(linSys[0].rows(), Lorder, W);
        weigths = ITHACAregularization::PCGLS(Theta, L, W, 
                linSys[1], CG_Nsteps);
    }
    else
    {
        Info << "Select a linear system solver in this list:" << endl
             << "fullPivLU, jacobiSvd, householderQr, ldlt, TSVD, Tikhonov, " <<  
             "conjugateGradient" << endl;
        exit(1);
    }
    return weigths;
}

void sequentialIHTP::set_valueFraction()
{
    fvMesh& mesh = _mesh();
    valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
    refGrad.resize(mesh.boundaryMesh()["coldSide"].size());
    Eigen::VectorXd faceCellDist =
        ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
    forAll (valueFraction, faceI)
    {
        valueFraction[faceI] = 1.0 / (1.0 + (thermalCond / HTC / faceCellDist(faceI)));
        refGrad[faceI] =  0.0;
    }
}

void sequentialIHTP::assignDirectBC(label timeI)
{
    fvMesh& mesh = _mesh();
    volScalarField& T = _T();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, RobinBC, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(T, patchI, - heatFlux[timeI] / thermalCond);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
}

void sequentialIHTP::solveT_ic(volScalarField _initialField)
{
    Info << "\nSolving FULL T_ic problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T_ic(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf;
    word outputFolder = "./ITHACAoutput/debugT_ic/";

    if(timeSampleI == 0)
    {
        ITHACAutilities::assignIF(T_ic, _initialField);
    }
    else
    {
        ITHACAutilities::assignIF(T_ic, Ttime[Ttime.size() - 1]);
    }

    T_ic_field.append(T_ic.clone());
    T_ic_time.resize(0);
    label timeI = 0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T_ic, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T_ic, patchI, homogeneousBC);
        }
    }

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T_ic) - fvm::laplacian(DT * diffusivity, T_ic)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T_ic);
        }

        T_ic_time.append(T_ic.clone());
        T_ic_field.append(T_ic.clone());
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T_ic_vector = fieldValueAtThermocouples(T_ic_time);
    T_ic_ready = 1;
    Info << "SolveT_ic ENDED\n" << endl;
}

void sequentialIHTP::getT_ic_modes()
{
    word outputFolder = "./ITHACAoutput/modes/";
    Info << "Computing " << NmodesT_ic << " T_ic modes" << endl;

    ITHACAPOD::getModes(T_ic_field, T_ic_modes, "T_ic",
                        0, 0, 0,
                        NmodesT_ic);
    PtrList<volScalarField> modes = T_ic_modes.toPtrList();
    forAll(modes, modeI)
    {
        ITHACAstream::exportSolution(modes[modeI], std::to_string(modeI + 1),
                outputFolder, "T_ic_modes");
    }
    Info << "T_ic modes COMPUTED\n" << endl;
}

void sequentialIHTP::projectT_ic()
{
    Info << "\n*****************************************************" << endl;
    Info << "Computing projection matrices" << endl;
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T_ic(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T_ic, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T_ic, patchI, homogeneousBC);
        }
    }

    Eigen::SparseMatrix<double> T_ic_implicitMatrix;
    Eigen::SparseMatrix<double> T_ic_explicitMatrix;

    ITHACAutilities::assignIF(T_ic, 1.0);
    fvScalarMatrix Teq(fvm::ddt(T_ic) - fvm::laplacian(DT * diffusivity, T_ic));
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Teq, T_ic_implicitMatrix, b);

    T_ic_modes.toEigen();
    T_ic_implicitMatrix_red = T_ic_modes.EigenModes[0].transpose() 
        * T_ic_implicitMatrix * T_ic_modes.EigenModes[0];
    T_ic_explicitMatrix_red = T_ic_modes.EigenModes[0].transpose() 
        * b.asDiagonal() * T_ic_modes.EigenModes[0];
    Info << "Projection matrices COMPUTED\n" << endl;
}

void sequentialIHTP::projectDirectOntoT_ic()
{
    /// Creation of the matrices to project direct solution at the last timestep
    /// onto the T_ic reduced space
    Info << "Computing direct problem projection matrices" << endl;
    int internalFieldSize = T_basis[0][0].internalField().size();
    Eigen::MatrixXd T_basis_Eigen(internalFieldSize, Nbasis);

    PtrList<volScalarField> T_basis_lastTime;
    M_Assert(T_basis[0].size() == NtimeStepsBetweenSamples, 
            "The basis for the direct problem have wrong dimention in time");

    forAll(T_basis, baseI)
    {
        volScalarField temp = T_basis[baseI][NtimeStepsBetweenSamples - 1];
        T_basis_lastTime.append(temp.clone());
    }

    T_basis_projectionMat = T_ic_modes.project(T_basis_lastTime); 
}

void sequentialIHTP::pointProjectionOffline()
{
    int Ncells = magicPoints.size();
    M_Assert(Ncells > 0, "Set the number of magic points");

    int lastTimestepID = NtimeStepsBetweenSamples - 1;

    pointsReconstructMatrix.resize(Ncells, NmodesT_ic);
    for(int cellI = 0; cellI < Ncells; cellI++)
    {
        pointsReconstructMatrix.row(cellI) = 
            T_ic_modes.EigenModes[0].row(magicPoints[cellI]);
    }
    pointT_basis_reconstructionMat = pointsReconstructMatrix * T_basis_projectionMat;
}

void sequentialIHTP::projectionErrorOffline()
{
    /// I compute the L2 norm of the T_basis and Tad perpendicular to the projection
    Info << "Computing the offline part for the projection error" << endl;

    //TODO You can consider all modes when computing the error and then 
    //choose the right ammount of modes based on the behaviour of the error
    int lastTimestepID = NtimeStepsBetweenSamples - 1;
    projectionErrorT_basis.resize(0);
    word outputFolder = "./ITHACAoutput/projectionError";


    forAll(T_basis, baseI)
    {
        volScalarField base = T_basis[baseI][lastTimestepID];
        volScalarField baseProj(base);
        T_ic_modes.projectSnapshot(base, baseProj, NmodesT_ic, "L2");
        volScalarField temp = base - baseProj;
        projectionErrorT_basis.append(temp.clone());
        ITHACAstream::exportSolution(projectionErrorT_basis[baseI], 
                std::to_string(baseI + 1),
                outputFolder, "projectionErrorT_basis");
    }
}

void sequentialIHTP::T_ic_offline(int NmagicPoints)
{
    getT_ic_modes();
    findMagicPoints(NmagicPoints);
    projectT_ic();
    projectDirectOntoT_ic();
    projectionErrorOffline();
    pointProjectionOffline();
}

void sequentialIHTP::solveDirect()
{
    M_Assert(offlineFlag, "Callable only during offline");
    M_Assert(diffusivity>1e-36, "Set the diffusivity value");

    restartOffline();
    volScalarField& T = _T();
    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();
    fv::options& fvOptions(_fvOptions());

    fvMesh& mesh = _mesh();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, RobinBC, refGrad, valueFraction);
        }
        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
        {
            ITHACAutilities::assignBC(T, patchI, - heatFlux[1] / thermalCond);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
    ITHACAutilities::assignIF(T, homogeneousBC);

    label timeI = 0;
    Ttime.resize(0);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT * diffusivity, T)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        Ttime.append(T.clone());

        runTime.printExecutionTime(Info);
        runTime.write();
        timeI++;
    }

    Info << "Direct computation ENDED" << endl;
    
}

void sequentialIHTP::readThermocouples()
{
    if (!thermocouplesRead)
    {
        Info << "Defining positions of thermocouples" << endl;
        word fileName = "./thermocouplesCellsID";

        if (ITHACAutilities::check_file(fileName + "_mat.txt"))
        {
            Info << "Reading thermocouples cells from file" << endl;
            Eigen::MatrixXi TCmatrix = 
                ITHACAstream::readMatrix(fileName + "_mat.txt").cast<int> ();
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
            ITHACAstream::exportSolution(thermocouplesField, "1", 
                    "./ITHACAoutput/thermocouplesField/", "thermocouplesField");
            Eigen::MatrixXi thermocouplesCellID_eigen = Foam2Eigen::List2EigenMatrix(
                        thermocouplesCellID);
            ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
                                       "eigen", "./");
        }

        thermocouplesRead = 1;
        samplingTime.resize(timeSamplesNum);
        forAll(samplingTime, timeI)
        {
            samplingTime[timeI] = timeSamplesT0 + timeI * timeSamplesDeltaT;
        }
        sampling2symulationTime();
	NtimeStepsBetweenSamples = timeSamplesDeltaT / deltaTime;
        residual.resize(thermocouplesNum * timeSamplesNum);
    }
    else
    {
        WarningInFunction << "readThermocouples function called twice." << endl;
        WarningInFunction << "I am not doing the second reading." << endl;
    }
}

scalar sequentialIHTP::fieldValueAtPoint(
    volScalarField& _field, vector _point)
{
    fvMesh& mesh = _mesh();
    volScalarField field
    (
        "T",
        _field
    );
    dictionary interpolationDict =
        mesh.solutionDict().subDict("interpolationSchemes");
    autoPtr<Foam::interpolation<scalar>> fieldInterp =
        Foam::interpolation<scalar>::New(interpolationDict, field);
    return fieldInterp->interpolate(_point, mesh.findCell(_point));
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    volScalarField field)
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
        fieldInt(tcI) = fieldInterp->interpolate(thermocouplesPos[tcI],
                        thermocouplesCellID[tcI]);
    }
    return fieldInt;
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    PtrList<volScalarField> fieldList, label fieldI)
{
    Eigen::VectorXd fieldInt = fieldValueAtThermocouples(fieldList[fieldI]);
    return fieldInt;
}

Eigen::VectorXd sequentialIHTP::fieldValueAtThermocouples(
    PtrList<volScalarField> fieldList)
{
    Eigen::VectorXd fieldInt;
    if ( fieldList.size() == Ntimes + 1 )
    {
        Info << "\n Sampling for ALL sampling times \n\n" << endl;
        fieldInt.resize(timeSamplesNum * thermocouplesNum);
        forAll(samplingSteps, sampleTimeI)
        {
            Eigen::VectorXd temp = 
                fieldValueAtThermocouples(fieldList, samplingSteps[sampleTimeI]);
            fieldInt.segment(sampleTimeI * thermocouplesNum, thermocouplesNum) = temp;
        }
    }
    else if ( fieldList.size() == NtimeStepsBetweenSamples )
    {
        Info << "\nSampling ONLY the last timestep\n\n" << endl;
        fieldInt =
            fieldValueAtThermocouples(fieldList, NtimeStepsBetweenSamples - 1);
    }
    else
    {
        Info << "The input fieldList of sequentialIHTP::fieldValueAtThermocouples " << 
            "can have size Ntimes + 1 (=" << Ntimes + 1 << ") or\n";
        Info << " NtimeStepsBetweenSamples  (=" <<  NtimeStepsBetweenSamples << 
            ") but has size " << fieldList.size() << endl;
        Info << "Exiting." << endl;
        exit(23);
    }
    Info << "\nSampling done \n" << endl;
    return fieldInt;
}


void sequentialIHTP::restart()
{
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(Times[1], 0);
    _simple.clear();
   _T.clear();

    Foam::fvMesh& mesh = _mesh();
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );

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
    

    Info << "Ready for new computation" << endl;
}

void sequentialIHTP::restartOffline()
{
    Info << "Setting endTime to offlineEndTime" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(offlineEndTime);
    Info << "Ready for new offline computation" << endl;
}

void sequentialIHTP::restartT_ic()
{
    Info << "Setting endTime to offlineEndTime" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(timeSamplesDeltaT);
    Info << "Ready for new T_ic computation" << endl;
}

void sequentialIHTP::sampling2symulationTime()
{
    scalar EPSILON = 2e-16;
    label deltaTimeQuotient = std::floor(timeSamplesDeltaT / deltaTime);
    //Info << "debug: deltaTimeQuotient= " <<deltaTimeQuotient<< endl;
    //Info << "debug: timeSamplesDeltaT= " <<timeSamplesDeltaT<< endl;
    //Info << "debug: deltaTime= " <<deltaTime<< endl;
    M_Assert(std::fabs(timeSamplesDeltaT / deltaTime - std::trunc(
                           timeSamplesDeltaT / deltaTime)) < EPSILON,
             "timeSamplesDeltaT should be a multiple of deltaTime");
    label n0 = (timeSamplesT0 - startTime) / deltaTime;
    M_Assert(n0 > 0, "First sampling step cannot be 0");
    //Info << "debug: n0 = " << n0 << endl;
    //Info << "debug: (timeSamplesT0 - startTime) / deltaTime = " << (timeSamplesT0 - startTime) / deltaTime << endl;
    M_Assert(std::fabs(n0 * deltaTime - timeSamplesT0) < EPSILON,
             "The first sampling time must coincide with a simulation timestep");
    scalar samplingEndTime = timeSamplesDeltaT * (timeSamplesNum - 1) + timeSamplesT0;
    //Info << "debug: samplingEndTime = " << samplingEndTime << endl;
    //Info << "debug: EndTime = " << endTime << endl;
    M_Assert(!(endTime + EPSILON < samplingEndTime
               && std::fabs(endTime - samplingEndTime) > EPSILON),
             "The samplingEndTime cannot be later than the symulation endTime");
    samplingSteps.resize(timeSamplesNum);
    forAll(samplingTime, sampleI)
    {
        samplingSteps[sampleI] = n0 + sampleI * deltaTimeQuotient;
    }
    Info << "debug: samplingSteps = " << samplingSteps << endl;
}

void sequentialIHTP::parameterizedHeatFlux_postProcess(
    List<Eigen::MatrixXd> linSys, Eigen::VectorXd weigths, word outputFolder, 
    label verbose)
{
    Foam::fvMesh& mesh = _mesh();
    Eigen::VectorXd Tcomp = fieldValueAtThermocouples(Ttime);
    //std::cout << "Tcomp = \n" << Tcomp.transpose() << std::endl;
    //std::cout << "TmeasShort = \n" << TmeasShort.transpose() << std::endl;
    J = 0.5 * Foam::sqrt((Tcomp - TmeasShort).dot(Tcomp - TmeasShort));
    Info << "J = " << J << endl;
    Jlist.conservativeResize(Jlist.size() + 1);
    Jlist(Jlist.size() - 1) = J;
    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", outputFolder);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Theta,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singVal = svd.singularValues();
    ITHACAstream::exportMatrix(singVal, "singularValues", "eigen", outputFolder);
    if (verbose)
    {
        // Printing outputs at screen
        std::cout << "Singular values of Theta.transpose() * Theta are " << std::endl;
        std::cout << svd.singularValues() << std::endl;
        std::cout << "weigths = " << std::endl;
        std::cout << weigths << std::endl;
        std::cout << "linSys[1] = " << std::endl;
        std::cout << linSys[1] << std::endl;
        std::cout << "Theta = " << std::endl;
        std::cout << Theta << std::endl;
        residual =  linSys[0] * weigths - linSys[1];
        //std::cout << "Residual  = " << std::endl;
        //std::cout << residual << std::endl;
        std::cout << "Residual 2-norm = " << std::endl;
        std::cout << residual.squaredNorm() << std::endl;
        std::cout << "T_ic_vector = " << std::endl;
        std::cout << T_ic_vector << std::endl;
        std::cout << "Tmeas = " << std::endl;
        std::cout << Tmeas << std::endl;
    }
}

void sequentialIHTP::findMagicPoints(int NmagicPoints)
{
    M_Assert(NmagicPoints > 0, "Set number of magic points");
    magicPoints.clear();
    M_Assert(NmagicPoints <= NmodesT_ic, 
            "Number of magic points bigger than number of modes");
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd r;
    Eigen::VectorXd rho(1);
    Eigen::MatrixXd MatrixModes = T_ic_modes.toEigen()[0];
    label ind_max, c1;
    double max = MatrixModes.cwiseAbs().col(0).maxCoeff(&ind_max, &c1);
    rho(0) = max;
    magicPoints.append(ind_max);
    Eigen::MatrixXd U = MatrixModes.col(0);
    Eigen::SparseMatrix<double> P;
    P.resize(MatrixModes.rows(), 1);
    P.insert(ind_max, 0) = 1;

    for (label i = 1; i < NmagicPoints; i++)
    {
        A = P.transpose() * U;
        b = P.transpose() * MatrixModes.col(i);
        c = A.fullPivLu().solve(b);
        r = MatrixModes.col(i) - U * c;
        max = r.cwiseAbs().maxCoeff(&ind_max, &c1);
        P.conservativeResize(MatrixModes.rows(), i + 1);
        P.insert(ind_max, i) = 1;
        U.conservativeResize(MatrixModes.rows(), i + 1);
        U.col(i) =  MatrixModes.col(i);
        rho.conservativeResize(i + 1);
        rho(i) = max;
        magicPoints.append(ind_max);
    }

    Info << "magicPoints:\n" << magicPoints << endl;
}

void sequentialIHTP::computeBasisCrossIntegral()
{
    Info << "Computing the cross L2 product of the space basis " << endl;
    basisCrossIntegralMatrix.resize(Nbasis, Nbasis);

    fvMesh& mesh = _mesh();
    M_Assert(heatFluxSpaceBasis.size() == Nbasis, 
            "heatFluxSpaceBasis are not properly set up");
    forAll(heatFluxSpaceBasis, baseI)
    {
        forAll(heatFluxSpaceBasis, baseJ)
        {
            basisCrossIntegralMatrix(baseI, baseJ) = 
                ITHACAutilities::L2productOnPatch(mesh, 
                        heatFluxSpaceBasis[baseI], heatFluxSpaceBasis[baseJ], "hotSide");
        }
    }

    Info << "computeBasisCrossIntegral END" << endl;
}
