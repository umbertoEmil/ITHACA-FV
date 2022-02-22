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
/// Source file of the reducedSequentialIHTP_incrementalPOD_constant class

#include "ReducedSequentialIHTP_incrementalPOD_constant.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSequentialIHTP_incrementalPOD_constant::reducedSequentialIHTP_incrementalPOD_constant()
{
}

reducedSequentialIHTP_incrementalPOD_constant::reducedSequentialIHTP_incrementalPOD_constant(int argc, 
        char* argv[])
    :
    sequentialIHTP_constant(argc, argv)
{
    
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedSequentialIHTP_incrementalPOD_constant::computeHeatFluxWeights(
        volScalarField _initialField, List<scalar> _initialWeights, word _outputFolder,  
        int _NmagicPoints, double _SVDtol, word _PODnorm)
{
    timeSampleI = 0;
    onlineCPUtime_avg = 0;
    onlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    offlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    bool recomputeLastStep = 0;
    NmagicPoints = _NmagicPoints;
    M_Assert(T_ic_projectionTol > 0.0, "Initialize T_ic_projectionTol");

    heatFluxWeights = _initialWeights;

    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    if (linSys_solver == "TSVD" || linSys_solver == "Tikhonov" || 
            linSys_solver == "conjugateGradient" || linSys_solver == "PCGLS")
    {
        linSys[0] = Theta;
    }
    else
    {
        linSys[0] = Theta.transpose() * Theta;
    }

    bool comingBackFlag = 0; // Set to 1 if I am coming back of one step
    while(timeSampleI < timeSamplesNum)
    {
        Info << "\n\n************************************************************" << endl;
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        bool doOnline = 0;

        if(T_ic_field.size() >= NmodesT_ic + 2)
        {
            if(computeT_ic_projectionError(Ttime[NtimeStepsBetweenSamples -1]) > 
                    T_ic_projectionTol)
            {
                if(previousWasReduced)
                {
                    Info << "\nComing back of one step\n" << endl;
                    timeSampleI--;
                    ITHACAutilities::assignIF(Ttime[NtimeStepsBetweenSamples -1], 
                            _initialField);
                }
            }
            else
            {
                ITHACAutilities::assignIF(_initialField, 
                        Ttime[NtimeStepsBetweenSamples -1]);
                doOnline = 1;
            }
        }

        if(doOnline)
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            Info << "Using REDUCED T_ic" << endl;
            onlineCountVec(timeSampleI) = 1;
            onlineWindowsVec.conservativeResize(onlineWindowsVec.size() + 1);
            onlineWindowsVec(onlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            bool useReducedInitialField = 0;
            solveT_ic_online(_initialField, useReducedInitialField);
            previousWasReduced = 1;
            Info << "T_ic computed" << endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " 
                << std::chrono::duration_cast<std::chrono::microseconds> 
                (end - begin).count() << "[microseconds]" << std::endl;
        }
        else
        {
            Info << "Using FULL T_ic" << endl;
            comingBackFlag = 0;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            offlineCountVec(timeSampleI) = 1;
            offlineWindowsVec.conservativeResize(offlineWindowsVec.size() + 1);
            offlineWindowsVec(offlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            previousWasReduced = 0;
            solveT_ic(_initialField);
            Info << "T_ic computed" << endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " 
                << std::chrono::duration_cast<std::chrono::microseconds> 
                (end - begin).count() << "[microseconds]" << std::endl;
            if(IPOD.rank == 0)
            {
                IPOD.tolleranceSVD = _SVDtol;
                IPOD.PODnorm = _PODnorm;
                M_Assert(_PODnorm == "L2" ||
                        _PODnorm == "Frobenius", 
                        "The PODnorm can be only L2 or Frobenius");
                IPOD.initialize(T_ic_field[0]);
            }
            else
            {
                label fistTimestep = T_ic_field.size() - NtimeStepsBetweenSamples - 1;
                IPOD.addSnapshot(T_ic_field[fistTimestep]);
            }
            IPOD.addSnapshot(T_ic_time);
            findMagicPoints();
            projectT_ic();

            projectDirectOntoT_ic();
            projectionErrorOffline();
            pointProjectionOffline();
        }

        TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, thermocouplesNum *
                NsamplesWindow);

        if (linSys_solver == "TSVD" || linSys_solver == "Tikhonov" || 
                linSys_solver == "conjugateGradient" || linSys_solver == "PCGLS")
        {
            linSys[1] = TmeasShort - T_ic_vector;
        }
        else
        {
            linSys[1] = Theta.transpose() * ( TmeasShort - T_ic_vector );
        }

        Eigen::VectorXd weigths = solveLinSys(linSys);

        forAll(heatFluxWeights, weightI)
        {
            heatFluxWeights[weightI] = weigths(weightI);
        }
        updateHeatFlux(heatFluxWeights);

        label verbose = 0;
        reconstructT(_outputFolder);
        parameterizedHeatFlux_postProcess(linSys, weigths, _outputFolder, verbose);
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms =
            std::chrono::duration<double, std::milli>(t_end-t_start).count();
        onlineCPUtime_avg += elapsed_time_ms;
        Info << "CPU time = " << elapsed_time_ms << " milliseconds" << endl << endl;

        timeSampleI++;
    }
    onlineCPUtime_avg = onlineCPUtime_avg / timeSamplesNum;
    Eigen::VectorXd onlineCPUtime_eig(1);
    onlineCPUtime_eig(0) = onlineCPUtime_avg;
    ITHACAstream::exportMatrix(
                    onlineCPUtime_eig, "onlineCPUtime_avg", "eigen", "./");

    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", _outputFolder);
    ITHACAstream::exportMatrix(onlineCountVec, "onlineCountVec", "eigen", _outputFolder);
    ITHACAstream::exportMatrix(onlineWindowsVec, "onlineWindowsVec", "eigen", 
            _outputFolder);
    ITHACAstream::exportMatrix(offlineCountVec, "offlineCountVec", "eigen", _outputFolder);
    ITHACAstream::exportMatrix(offlineWindowsVec, "offlineWindowsVec", "eigen", 
            _outputFolder);
    Info << "End" << endl;
    Info << endl;
}

void reducedSequentialIHTP_incrementalPOD_constant::solveT_ic_online(
        volScalarField _initialField, bool _useReducedInitialField)
{
    Info << "\nSolving REDUCED T_ic problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T_ic(_T);
    Foam::Time& runTime = _runTime();

    if(_useReducedInitialField)
    {
        Eigen::VectorXd heatFluxWeights_Eig = 
            Foam2Eigen::List2EigenMatrix(heatFluxWeights);

        if(T_ic_red.size() == 0 || !previousWasReduced)
        {
            T_ic_red = IPOD.project(_initialField);
        }
        else
        {
            T_ic_red = T_basis_projectionMat * heatFluxWeights_Eig + T_ic_red;
        }
    }
    else
    {
        T_ic_red = IPOD.project(_initialField);
    }

    T_ic_time.resize(0);
    label timeI = 0;
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        /// Reduced
        T_ic_red = T_ic_implicitMatrix_red.fullPivLu().solve(
                T_ic_explicitMatrix_red * T_ic_red + T_ic_source_red); 
        
        T_ic = IPOD.reconstruct(T_ic, T_ic_red, "T");
        T_ic_time.append(T_ic.clone());

        
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T_ic_vector = fieldValueAtThermocouples(T_ic_time);
    T_ic_ready = 1;

    Info << "SolveT_ic ENDED\n" << endl;
}

double reducedSequentialIHTP_incrementalPOD_constant::computeT_ic_projectionError(
        volScalarField& T_ic_in)
{
    double error = 0;

    // Pointwise
    Eigen::VectorXd TatPoints(magicPoints.size());
    Eigen::VectorXd TperpAtPoints = TatPoints;
    Eigen::VectorXd TprojAtPoints = TatPoints;
    Eigen::VectorXd magicPointsVolume = TatPoints;
    fvMesh& mesh = _mesh();

    if(previousWasReduced)
    {
        Info <<"Previous step was reduced" << endl;
        TprojAtPoints = pointT_basis_reconstructionMat * 
            Foam2Eigen::List2EigenMatrix(heatFluxWeights)
            + pointsReconstructMatrix * T_ic_red; 
        forAll(magicPoints, cellI)
        {
            TatPoints(cellI) = T_ic_in.internalField()[magicPoints[cellI]];
        }
    }
    else
    {
        Info << "Previous step was full" << endl;
        volScalarField T_ic_proj(IPOD.projectSnapshot(T_ic_in));                                
        forAll(magicPoints, cellI)
        {
            TprojAtPoints(cellI) = T_ic_proj.internalField()[magicPoints[cellI]];
            TatPoints(cellI) = T_ic_in.internalField()[magicPoints[cellI]];
        }
    }


    error = 0;
    TperpAtPoints = TatPoints - TprojAtPoints;
    forAll(magicPoints, cellI)
    {
        double tmp = TperpAtPoints(cellI) / TatPoints(cellI);
        error += tmp * tmp * mesh.V()[magicPoints[cellI]];
    }
    std::cout << "Projection error Points = " << error << std::endl;

    return error;
}

void reducedSequentialIHTP_incrementalPOD_constant::findMagicPoints()
{
    magicPoints.clear();
    NmagicPoints = IPOD.rank;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd r;
    Eigen::VectorXd rho(1);
    Eigen::MatrixXd MatrixModes = IPOD.toEigen()[0];
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

void reducedSequentialIHTP_incrementalPOD_constant::projectT_ic()
{
    Info << "\n*****************************************************" << endl;
    Info << "Computing projection matrices" << endl;
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T0(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T0, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T0, patchI, homogeneousBC);
        }
    }

    Eigen::SparseMatrix<double> T_ic_implicitMatrix;
    Eigen::SparseMatrix<double> T_ic_explicitMatrix;

    ITHACAutilities::assignIF(T0, 1.0);
    fvScalarMatrix Teq(fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0));
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Teq, T_ic_implicitMatrix, b);

    Info << "debug : b.sum() = " << b.sum() << endl;

    T_ic_implicitMatrix_red = IPOD.EigenModes[0].transpose()
        * T_ic_implicitMatrix * IPOD.EigenModes[0];
    T_ic_explicitMatrix_red = IPOD.EigenModes[0].transpose()
        * b.asDiagonal() * IPOD.EigenModes[0];
    Info << "debug : T_ic_explicitMatrix_red.sum() = " << T_ic_explicitMatrix_red.sum() << endl;

    RobinBC = Tf;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T0, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T0, patchI, homogeneousBC);
        }
    }

    ITHACAutilities::assignIF(T0, 1.0);
    fvScalarMatrix Teq2(fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0));
    Eigen::VectorXd b2;
    Foam2Eigen::fvMatrix2Eigen(Teq2, T_ic_implicitMatrix, b2);

    Info << "debug : b2.sum() = " << b2.sum() << endl;
    Eigen::VectorXd bDiff = b2 - b;
    Info << "debug : bDiff.sum() = " << bDiff.sum() << endl;

    T_ic_source_red = IPOD.EigenModes[0].transpose()
        * bDiff;
    Info << "debug : T_ic_source_red.sum() = " << T_ic_source_red.sum() << endl;
    Info << "Projection matrices COMPUTED\n" << endl;
}


void reducedSequentialIHTP_incrementalPOD_constant::projectDirectOntoT_ic()
{
    /// Creation of the matrices to project direct solution at the last timestep
    /// onto the T0 reduced space
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

    T_basis_projectionMat = IPOD.project(T_basis_lastTime);
}

void reducedSequentialIHTP_incrementalPOD_constant::projectionErrorOffline()
{
    /// I compute the L2 norm of the T_basis perpendicular to the projection
    Info << "Computing the offline part for the projection error" << endl;

    //TODO You can consider all modes when computing the error and then
    //choose the right ammount of modes based on the behaviour of the error
    int lastTimestepID = NtimeStepsBetweenSamples - 1;
    projectionErrorT_basis.resize(0);
    word outputFolder = "./ITHACAoutput/projectionError";

    forAll(T_basis, baseI)
    {
        volScalarField base = T_basis[baseI][lastTimestepID];
        volScalarField baseProj(IPOD.projectSnapshot(base));
        volScalarField temp = base - baseProj;
        projectionErrorT_basis.append(temp.clone());
        ITHACAstream::exportSolution(projectionErrorT_basis[baseI], 
                std::to_string(baseI + 1), outputFolder, "projectionErrorT_basis");
    }
    Info << "projectionErrorOffline END" << endl;
}

void reducedSequentialIHTP_incrementalPOD_constant::pointProjectionOffline()
{
    Info << "pointProjectionOffline START" << endl;
    int Ncells = magicPoints.size();
    M_Assert(Ncells > 0, "Set the number of magic points");

    int lastTimestepID = NtimeStepsBetweenSamples - 1;

    pointsReconstructMatrix.resize(Ncells, IPOD.rank);
    for(int cellI = 0; cellI < Ncells; cellI++)
    {
        pointsReconstructMatrix.row(cellI) =
            IPOD.EigenModes[0].row(magicPoints[cellI]);
    }
    pointT_basis_reconstructionMat = pointsReconstructMatrix * T_basis_projectionMat;
    Info << "pointProjectionOffline END" << endl;
}




