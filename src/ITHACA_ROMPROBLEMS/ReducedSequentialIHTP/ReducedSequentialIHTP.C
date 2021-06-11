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
/// Source file of the reducedSequentialIHTP class

#include "ReducedSequentialIHTP.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSequentialIHTP::reducedSequentialIHTP()
{
}

reducedSequentialIHTP::reducedSequentialIHTP(int argc, char* argv[])
    :
    sequentialIHTP(argc, argv)
{
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

void reducedSequentialIHTP::parameterizedBC(word outputFolder, volScalarField initialField,
        int NmagicPoints, word linSys_solver, label TSVD_filter)
{
    Info << endl << "Using quasilinearity of direct problem AND reduced T0" << endl;
    timeSampleI = 0;
    onlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    offlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
    bool recomputeLastStep = 0;
    M_Assert(T0projectionTol > 0.0, "Initialize T0projectionTol");

    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        Info << "debug: T0field.size() = " << T0field.size() << endl; 
        volScalarField oldInitialField = initialField;
        if(timeSampleI > 0)
        {
            /// Assign the new initialField
            reconstrucT("./ITHACAoutput/debugReconstrucT/");
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
                << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count() 
                << "[microseconds]" << std::endl;
        }
        else
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            offlineCountVec(timeSampleI) = 1;
            offlineWindowsVec.conservativeResize(offlineWindowsVec.size() + 1);
            offlineWindowsVec(offlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            previousWasReduced = 0;
            Info << "Using FULL T0" << endl;
            solveT0(initialField);
            Info << "T0 computed" << endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Time difference = " 
                << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count() 
                << "[microseconds]" << std::endl;
            if(T0field.size() >= NmodesT0 + 2)
            {
                T0offline(NmagicPoints);
            }
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
            weigths = ITHACAregularization::TSVD(linSys[0], linSys[1], TSVD_filter);
        }
        else
        {
            Info << "Select a linear system solver in this list:" << endl
                 << "fullPivLU, jacobiSvd, householderQr, ldlt" << endl;
            exit(1);
        }
        gWeightsOld = gWeights;
        gWeights.resize(weigths.size());
        forAll(gWeights, weightI)
        {
            gWeights[weightI] = weigths(weightI);
        }
        Info << "Weights = \n" << gWeights << endl;
        update_gParametrized(gWeights);
        label verbose = 0;
        parameterizedBC_postProcess(linSys, weigths, outputFolder, verbose);
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

void reducedSequentialIHTP::solveT0online(volScalarField initialField, 
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
        Eigen::VectorXd gWeights_Eig = Foam2Eigen::List2EigenMatrix(gWeights);

        if(T0red.size() == 0 || !previousWasReduced)
        {
            T0red = T0modes.project(initialField);
        }
        else
        {
            T0red = Tbasis_projectionMat * gWeights_Eig - Tad_projected + T0red;
        }
    }
    else
    {
        ITHACAutilities::assignIF(T0, initialField);
        T0red = T0modes.project(T0);
    }

    T0_time.resize(0);
    label timeI = 0;
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        /// Reduced
        T0red = T0implicitMatrix_red.fullPivLu().solve(T0explicitMatrix_red * T0red); 
        
        T0 = T0modes.reconstruct(T0, T0red, "T");
        T0_time.append(T0.clone());

        
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T0_vector = fieldValueAtThermocouples(T0_time);
    Info << "SolveT0 ENDED\n" << endl;
}

double reducedSequentialIHTP::T0projectionError(volScalarField& T0in)
{
    double error = 0;
    //int lastTimestepID = NtimeStepsBetweenSamples - 1;
    //volScalarField errorField(_T());
    //ITHACAutilities::assignIF(errorField, 0.0);
    //forAll(Tbasis, baseI)                                                              
    //{                                                                                  
    //    volScalarField temp (projectionErrorTbasis[baseI] + projectionErrorTad[0]);
    //    error += std::abs(gWeights[baseI]) * ITHACAutilities::L2Norm(temp);
    //}                                                                                  
    //if(previousWasReduced)
    //{
    //    error += ITHACAutilities::L2Norm(projectionErrorTad[0]);
    //}
    //else
    //{
    //    volScalarField T0Proj(T0_time[lastTimestepID]);                                 
    //    T0modes.projectSnapshot(T0_time[lastTimestepID], T0Proj, NmodesT0, "L2");
    //    volScalarField T0Perp = T0_time[lastTimestepID] - T0Proj;
    //    error += ITHACAutilities::L2Norm(projectionErrorTad[0]) + 
    //        ITHACAutilities::L2Norm(T0Perp);               
    //}
    //error = std::abs(error);

    //Info << "Error Norm Reduced = " << error << endl;
    //Info << "Projection Reduced = " << error / ITHACAutilities::L2Norm(T0in) << endl;

    //volScalarField T0inProjected(T0in);
    //T0modes.projectSnapshot(T0in, T0inProjected, NmodesT0, "L2");
    //volScalarField T0inPerp(T0in - T0inProjected);
    //error = ITHACAutilities::L2Norm(T0inPerp) / ITHACAutilities::L2Norm(T0in);
    //Info << "Error Norm = " << ITHACAutilities::L2Norm(T0inPerp) << endl;
    //Info << "Projection error = " << error << endl;

    // Pointwise
    Eigen::VectorXd TatPoints(magicPoints.size());
    Eigen::VectorXd TperpAtPoints = TatPoints;
    Eigen::VectorXd TprojAtPoints = TatPoints;
    Eigen::VectorXd magicPointsVolume = TatPoints;
    fvMesh& mesh = _mesh();

    //double error_num = 0;
    //double error_den = 0;
    //error = 0;
    //forAll(magicPoints, cellI)
    //{
    //    TatPoints(cellI) = T0in.internalField()[magicPoints[cellI]];
    //    TprojAtPoints(cellI) = T0inProjected.internalField()[magicPoints[cellI]];
    //    TperpAtPoints(cellI) = T0inPerp.internalField()[magicPoints[cellI]];
    //    magicPointsVolume(cellI) = mesh.V()[magicPoints[cellI]];
    //    error_num += TperpAtPoints(cellI) * TperpAtPoints(cellI) * 
    //        magicPointsVolume(cellI);
    //    error_den += TatPoints(cellI) * TatPoints(cellI) * magicPointsVolume(cellI);
    //    double tmp = TperpAtPoints(cellI) / TatPoints(cellI);
    //    error += tmp * tmp * magicPointsVolume(cellI);
    //}


    if(previousWasReduced)
    {
        Info <<"Previous step was reduced" << endl;
        TprojAtPoints = pointTbasis_reconstructionMat * 
            Foam2Eigen::List2EigenMatrix(gWeights) - pointTad_reconstructed 
            + pointsReconstructMatrix * T0red; 
        forAll(magicPoints, cellI)
        {
            TatPoints(cellI) = T0in.internalField()[magicPoints[cellI]];
        }
    }
    else
    {
        Info << "Previous step was full" << endl;
        volScalarField T0Proj(T0in);                                
        T0modes.projectSnapshot(T0in, T0Proj, NmodesT0, "L2");
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

volScalarField reducedSequentialIHTP::reconstrucT_lastTime()
{                                                                                          
    Info << "Reconstructing last timestesp field T" << endl;
    Ttime.resize(0);                                                                       
    restart();

    int timeI = NtimeStepsBetweenSamples - 1;
    volScalarField Tout(_T);
    ITHACAutilities::assignIF(Tout, homogeneousBC);                                       
    forAll(Tbasis, baseI)                                                              
    {                                                                                  
        Tout += gWeights[baseI] * (Tbasis[baseI][timeI] + Tad_time[timeI]);
    }                                                                                  
    Tout += - Tad_time[timeI] + T0_time[timeI];
    return Tout;
}   
