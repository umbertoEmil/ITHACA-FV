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
/// Source file of the reducedSequentialIHTP_incrementalPOD class

#include "ReducedSequentialIHTP_incrementalPOD.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSequentialIHTP_incrementalPOD::reducedSequentialIHTP_incrementalPOD()
{
}

reducedSequentialIHTP_incrementalPOD::reducedSequentialIHTP_incrementalPOD(int argc, 
        char* argv[])
    :
    sequentialIHTP(argc, argv)
{
    
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedSequentialIHTP_incrementalPOD::parameterizedBC(word outputFolder, 
        volScalarField initialField, int _NmagicPoints, double _SVDtol, word _PODnorm)
{
    Info << endl << "Using quasilinearity of direct problem AND reduced T0" << endl;
    Info << "Using " << linSys_solver << " to solve the linear system" << endl;
    timeSampleI = 0;
    onlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    offlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    bool recomputeLastStep = 0;
    NmagicPoints = _NmagicPoints;
    M_Assert(T0projectionTol > 0.0, "Initialize T0projectionTol");

    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        volScalarField oldInitialField = initialField;
        if(timeSampleI > 0)
        {
            /// Assign the new initialField
            reconstrucT(initialField, "./ITHACAoutput/debugReconstrucT/");
            ITHACAutilities::assignIF(initialField, Ttime[NtimeStepsBetweenSamples -1]);
        }

        bool doOnline = 0;
        if(T0field.size() >= NmodesT0 + 2)
        {
            if(T0projectionError(initialField) > T0projectionTol)
            {
                if(previousWasReduced)
                {
                    Info << "\ndebug: Coming back of one step\n" << endl;
                    timeSampleI--;
                    ITHACAutilities::assignIF(initialField, oldInitialField);
                }
            }
            else
            {
                doOnline = 1;
            }
        }

        if(doOnline)
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            Info << "Using REDUCED T0" << endl;
            onlineCountVec(timeSampleI) = 1;
            onlineWindowsVec.conservativeResize(onlineWindowsVec.size() + 1);
            onlineWindowsVec(onlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            bool useReducedInitialField = 1;
            solveT0online(initialField, useReducedInitialField);
            previousWasReduced = 1;
            Info << "T0 computed" << endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " 
                << std::chrono::duration_cast<std::chrono::microseconds> 
                (end - begin).count() << "[microseconds]" << std::endl;
        }
        else
        {
            Info << "Using FULL T0" << endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            offlineCountVec(timeSampleI) = 1;
            offlineWindowsVec.conservativeResize(offlineWindowsVec.size() + 1);
            offlineWindowsVec(offlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            previousWasReduced = 0;
            solveT0(initialField);
            Info << "T0 computed" << endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " 
                << std::chrono::duration_cast<std::chrono::microseconds> 
                (end - begin).count() << "[microseconds]" << std::endl;
            if(IPOD.rank == 0)
            {
                IPOD.tolleranceSVD = _SVDtol;
                IPOD.PODnorm = _PODnorm;
                M_Assert(_PODnorm == "L2" ||
                        _PODnorm == "Frobenius", "The PODnorm can be only L2 or Frobenius");
                IPOD.initialize(T0field[0]);
            }
            else
            {
                label fistTimestep = T0field.size() - NtimeStepsBetweenSamples - 1;
                IPOD.addSnapshot(T0field[fistTimestep]);
            }
            IPOD.addSnapshot(T0_time);
            findMagicPoints();
            projectT0();

            projectDirectOntoT0();
            projectionErrorOffline();
            pointProjectionOffline();
        }
        List<Eigen::MatrixXd> linSys;
        linSys.resize(2);
                                                                                           
        TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, 
                thermocouplesNum * NsamplesWindow);
        linSys[0] = Theta.transpose() * Theta;                                             
        linSys[1] = Theta.transpose() * (TmeasShort + addSol - T0_vector);                 
        Eigen::VectorXd weigths;                                                           
        
        if (linSys_solver == "fullPivLU")                                                  
        {
            weigths = linSys[0].fullPivLu().solve(linSys[1]);                              
        }
        else if (linSys_solver == "jacobiSvd")
        {                                                                                  
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                                  Eigen::ComputeThinU | Eigen::ComputeThinV
                                                  );           
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
            weigths = ITHACAregularization::TSVD(linSys[0], linSys[1], TSVD_filter);
        }
        else if (linSys_solver == "Tikhonov")
        {
            weigths = ITHACAregularization::Tikhonov(linSys[0], linSys[1], Tikhonov_filter);
        }
        else
        {
            Info << "Select a linear system solver in this list:" << endl
                 << "fullPivLU, jacobiSvd, householderQr, ldlt, TSVD, Tikhonov" << endl;
            exit(1);
        }

        heatFluxWeightsOld = heatFluxWeights;
        heatFluxWeights.resize(weigths.size());
        forAll(heatFluxWeights, weightI)
        {
            heatFluxWeights[weightI] = weigths(weightI);
        }
        Info << "Weights = \n" << heatFluxWeights << endl;
        updateHeatFlux(heatFluxWeights);
        label verbose = 0;
        parameterizedBC_postProcess(linSys, weigths, initialField, outputFolder, verbose);
        timeSampleI++;
    }
    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", outputFolder);
    ITHACAstream::exportMatrix(onlineCountVec, "onlineCountVec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(onlineWindowsVec, "onlineWindowsVec", "eigen", 
            outputFolder);
    ITHACAstream::exportMatrix(offlineCountVec, "offlineCountVec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(offlineWindowsVec, "offlineWindowsVec", "eigen", 
            outputFolder);
    Info << "End" << endl;
    Info << endl;
}

void reducedSequentialIHTP_incrementalPOD::solveT0online(volScalarField initialField, 
        bool useReducedInitialField)
{
    Info << "\nSolving REDUCED T0 problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T0(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;
    word outputFolder = "./ITHACAoutput/debugT0/";
    Eigen::VectorXd T0red_prova;

    if(useReducedInitialField)
    {
        Eigen::VectorXd heatFluxWeights_Eig = Foam2Eigen::List2EigenMatrix(heatFluxWeights);

        if(T0red.size() == 0 || !previousWasReduced)
        {
            T0red = IPOD.project(initialField);
        }
        else
        {
            T0red = Tbasis_projectionMat * heatFluxWeights_Eig - Tad_projected + T0red;
        }
    }
    else
    {
        ITHACAutilities::assignIF(T0, initialField);
        T0red = IPOD.project(T0);
    }

    T0_time.resize(0);
    label timeI = 0;
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        /// Reduced
        T0red = T0implicitMatrix_red.fullPivLu().solve(T0explicitMatrix_red * T0red); 
        
        T0 = IPOD.reconstruct(T0, T0red, "T");
        T0_time.append(T0.clone());

        
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T0_vector = fieldValueAtThermocouples(T0_time);
    Info << "SolveT0 ENDED\n" << endl;
}

double reducedSequentialIHTP_incrementalPOD::T0projectionError(volScalarField& T0in)
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
        TprojAtPoints = pointTbasis_reconstructionMat * 
            Foam2Eigen::List2EigenMatrix(heatFluxWeights) - pointTad_reconstructed 
            + pointsReconstructMatrix * T0red; 
        forAll(magicPoints, cellI)
        {
            TatPoints(cellI) = T0in.internalField()[magicPoints[cellI]];
        }
    }
    else
    {
        Info << "Previous step was full" << endl;
        volScalarField T0Proj(IPOD.projectSnapshot(T0in));                                
        forAll(magicPoints, cellI)
        {
            TprojAtPoints(cellI) = T0Proj.internalField()[magicPoints[cellI]];
            TatPoints(cellI) = T0in.internalField()[magicPoints[cellI]];
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

volScalarField reducedSequentialIHTP_incrementalPOD::reconstrucT_lastTime()
{                                                                                          
    Info << "Reconstructing last timestesp field T" << endl;
    Ttime.resize(0);                                                                       
    restart();

    int timeI = NtimeStepsBetweenSamples - 1;
    volScalarField Tout(_T);
    ITHACAutilities::assignIF(Tout, homogeneousBC);                                       
    forAll(Tbasis, baseI)                                                              
    {                                                                                  
        Tout += heatFluxWeights[baseI] * (Tbasis[baseI][timeI] + Tad_time[timeI]);
    }                                                                                  
    Tout += - Tad_time[timeI] + T0_time[timeI];
    return Tout;
}   

void reducedSequentialIHTP_incrementalPOD::findMagicPoints()
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

void reducedSequentialIHTP_incrementalPOD::projectT0()
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

    Eigen::SparseMatrix<double> T0implicitMatrix;
    Eigen::SparseMatrix<double> T0explicitMatrix;

    ITHACAutilities::assignIF(T0, 1.0);
    fvScalarMatrix Teq(fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0));
    Eigen::VectorXd b;
    Foam2Eigen::fvMatrix2Eigen(Teq, T0implicitMatrix, b);

    T0implicitMatrix_red = IPOD.EigenModes[0].transpose()
        * T0implicitMatrix * IPOD.EigenModes[0];
    T0explicitMatrix_red = IPOD.EigenModes[0].transpose()
        * b.asDiagonal() * IPOD.EigenModes[0];
    Info << "Projection matrices COMPUTED\n" << endl;
}


void reducedSequentialIHTP_incrementalPOD::projectDirectOntoT0()
{
    /// Creation of the matrices to project direct solution at the last timestep
    /// onto the T0 reduced space
    Info << "Computing direct problem projection matrices" << endl;
    int internalFieldSize = Tbasis[0][0].internalField().size();
    Eigen::MatrixXd Tbasis_Eigen(internalFieldSize, Nbasis);
    Eigen::VectorXd Tad_Eigen =
        Foam2Eigen::field2Eigen(Tad_time[NtimeStepsBetweenSamples - 1]);

    PtrList<volScalarField> Tbasis_lastTime;
    M_Assert(Tbasis[0].size() == NtimeStepsBetweenSamples,
            "The basis for the direct problem have wrong dimention in time");

    forAll(Tbasis, baseI)
    {
        volScalarField temp = Tbasis[baseI][NtimeStepsBetweenSamples - 1]
            + Tad_time[NtimeStepsBetweenSamples - 1];
        Tbasis_lastTime.append(temp.clone());
    }

    Tbasis_projectionMat = IPOD.project(Tbasis_lastTime);
    Tad_projected = IPOD.project(Tad_time[NtimeStepsBetweenSamples - 1]);
}

void reducedSequentialIHTP_incrementalPOD::projectionErrorOffline()
{
    /// I compute the L2 norm of the Tbasis and Tad perpendicular to the projection
    Info << "Computing the offline part for the projection error" << endl;

    //TODO You can consider all modes when computing the error and then
    //choose the right ammount of modes based on the behaviour of the error
    int lastTimestepID = NtimeStepsBetweenSamples - 1;
    projectionErrorTbasis.resize(0);
    word outputFolder = "./ITHACAoutput/projectionError";

    forAll(Tbasis, baseI)
    {
        volScalarField base = Tbasis[baseI][lastTimestepID];
        volScalarField baseProj(IPOD.projectSnapshot(base));
        volScalarField temp = base - baseProj;
        projectionErrorTbasis.append(temp.clone());
        ITHACAstream::exportSolution(projectionErrorTbasis[baseI], 
                std::to_string(baseI + 1), outputFolder, "projectionErrorTbasis");
    }
    projectionErrorTad.resize(0);
    volScalarField TadProj(IPOD.projectSnapshot(Tad_time[lastTimestepID]));
    volScalarField temp(Tad_time[lastTimestepID] - TadProj);
    projectionErrorTad.append(temp.clone());
    ITHACAstream::exportSolution(projectionErrorTad[0], std::to_string(1),
            outputFolder, "projectionErrorTad");
}

void reducedSequentialIHTP_incrementalPOD::pointProjectionOffline()
{
    int Ncells = magicPoints.size();
    M_Assert(Ncells > 0, "Set the number of magic points");

    int lastTimestepID = NtimeStepsBetweenSamples - 1;

    pointsReconstructMatrix.resize(Ncells, IPOD.rank);
    for(int cellI = 0; cellI < Ncells; cellI++)
    {
        pointsReconstructMatrix.row(cellI) =
            IPOD.EigenModes[0].row(magicPoints[cellI]);
    }
    pointTbasis_reconstructionMat = pointsReconstructMatrix * Tbasis_projectionMat;
    pointTad_reconstructed = pointsReconstructMatrix * Tad_projected;
}




