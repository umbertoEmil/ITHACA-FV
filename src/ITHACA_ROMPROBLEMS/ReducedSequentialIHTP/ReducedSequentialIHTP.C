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
        Eigen::VectorXi errorCells, word linSys_solver, label TSVD_filter)
{
    Info << endl << "Using quasilinearity of direct problem AND reduced T0" << endl;
    timeSampleI = 0;
    onlineCountVec = Eigen::VectorXi::Zero(samplingTime.size());
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
            if(T0projectionError(initialField, errorCells) > T0projectionTol)
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
            Info << "Using REDUCED T0" << endl;
            onlineCountVec(timeSampleI) = 1;
            onlineWindowsVec.conservativeResize(onlineWindowsVec.size() + 1);
            onlineWindowsVec(onlineWindowsVec.size() - 1) = samplingTime[timeSampleI];
            bool useReducedInitialField = 1;
            solveT0online(initialField, 0);//, useReducedInitialField);
            previousWasReduced = 1;
            Info << "T0 computed" << endl;
        }
        else
        {
            previousWasReduced = 0;
            Info << "Using FULL T0" << endl;
            solveT0(initialField);
            if(T0field.size() >= NmodesT0 + 2)
            {
                T0offline(errorCells);
            }
        }
        List<Eigen::MatrixXd> linSys;
        linSys.resize(2);
                                                                                           
        TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, thermocouplesNum * NsamplesWindow);
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

    if(useReducedInitialField)
    {
        Eigen::VectorXd gWeights_Eig = Foam2Eigen::List2EigenMatrix(gWeights);

        if(T0red.size() == 0)
        {
            T0red = T0modes.project(initialField);
        }
        else
        {
            T0red = Tbasis_projectionMat * gWeights_Eig - Tad_projected + T0red;
        }
        //TODO projection error
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
        //std::cout << "debug : T0red = \n" << T0red << std::endl;
        
        T0 = T0modes.reconstruct(T0, T0red, "T");

        ITHACAstream::exportSolution(T0, std::to_string(timeSteps[samplingSteps[timeSampleI] - NtimeStepsBetweenSamples + timeI]),
                                     outputFolder,
                                     "T0reduced");
        T0_time.append(T0.clone());
        
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T0_vector = fieldValueAtThermocouples(T0_time);
    Info << "SolveT0 ENDED\n" << endl;
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

double reducedSequentialIHTP::T0projectionError(volScalarField& T0in, 
        Eigen::VectorXi errorCells)
{
    double error = 0;
    int lastTimestepID = NtimeStepsBetweenSamples - 1;
    volScalarField errorField(_T());
    ITHACAutilities::assignIF(errorField, 0.0);
    forAll(Tbasis, baseI)                                                              
    {                                                                                  
        //errorField += gWeights[baseI] * (projectionErrorTbasis[baseI] + projectionErrorTad[0]);
        volScalarField temp (projectionErrorTbasis[baseI] + projectionErrorTad[0]);
        error += std::abs(gWeights[baseI]) * ITHACAutilities::L2Norm(temp);
    }                                                                                  
    if(previousWasReduced)
    {
        error += ITHACAutilities::L2Norm(projectionErrorTad[0]);
    }
    else
    {
        volScalarField T0Proj(T0_time[lastTimestepID]);                                 
        T0modes.projectSnapshot(T0_time[lastTimestepID], T0Proj, NmodesT0, "L2");
        volScalarField T0Perp = T0_time[lastTimestepID] - T0Proj;
        error += ITHACAutilities::L2Norm(projectionErrorTad[0]) + ITHACAutilities::L2Norm(T0Perp);               
    }
    error = std::abs(error);

    Info << "Error Norm Reduced = " << error << endl;
    Info << "Projection Reduced = " << error / ITHACAutilities::L2Norm(T0in) << endl;
    volScalarField T0inProjected(T0in);
    T0modes.projectSnapshot(T0in, T0inProjected, NmodesT0, "L2");
    volScalarField T0inPerp(T0in - T0inProjected);
    error = ITHACAutilities::L2Norm(T0inPerp) / ITHACAutilities::L2Norm(T0in);
    Info << "Error Norm = " << ITHACAutilities::L2Norm(T0inPerp) << endl;
    Info << "Projection error = " << error << endl;


    // Pointwise
    Eigen::VectorXd TatPonts(errorCells.size());
    Eigen::VectorXd TperpAtPoint = TatPonts;
    Eigen::VectorXd TprojAtPoints = TatPonts;

    for(int cellI = 0; cellI < errorCells.size(); cellI++)
    {
        TatPonts(cellI) = T0in.internalField()[errorCells(cellI)];
    }

    if(previousWasReduced)
    {
        Info <<"Previous step was reduced" << endl;
        Info << "pointsProjectionMatrix = " << pointsProjectionMatrix.rows() << ", " << pointsProjectionMatrix.cols() << endl;
        Info << "T0red.size() = " << T0red.size() << endl;
        TprojAtPoints = pointsProjectionMatrix * T0red; 
    }
    else
    {
        Info <<"Previous step was full" << endl;
        volScalarField T0Proj(T0_time[lastTimestepID]);                                
        T0modes.projectSnapshot(T0_time[lastTimestepID], T0Proj, NmodesT0, "L2");
        for(int cellI = 0; cellI < errorCells.size(); cellI++)
        {
            TprojAtPoints(cellI) = T0Proj.internalField()[errorCells(cellI)];
        }
    }

    TperpAtPoint = TatPonts - TprojAtPoints;
    std::cout << "Projection error Points = " << TperpAtPoint.norm() / TatPonts.norm() 
        << std::endl;


    return error;
}
