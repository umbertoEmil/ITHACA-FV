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
/// Source file of the sequentialIHTP_constant class.


#include "sequentialIHTP_constant.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructors
sequentialIHTP_constant::sequentialIHTP_constant() {}

sequentialIHTP_constant::sequentialIHTP_constant(int argc, char* argv[])
    : 
    sequentialIHTP(argc, argv)
{
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

void sequentialIHTP_constant::updateTimeHeatFlux(List<scalar> weights)
{
    M_Assert(weights.size() == Nbasis,
             "weigths size different from basis functions size");
    forAll (weights, weightI)
    {
        for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
        {
            heatFluxTimeBasis[timeI][weightI] = weights[weightI];
        }
    }
}

void sequentialIHTP_constant::offlinePhase(bool force)
{
    fvMesh& mesh = _mesh();
    T_basis.resize(0);
    T_ic_field.resize(0);
    Theta.resize(thermocouplesNum, Nbasis);
    M_Assert(diffusivity > 0.0, "Call setDiffusivity to set up the diffusivity");

    offlineTimestepsSize = NtimeStepsBetweenSamples;
    offlineEndTime = NtimeStepsBetweenSamples * deltaTime;

    if (ITHACAutilities::check_file(folderOffline + "/Theta_mat.txt") && force == 0)
    {
        Info << "\nOffline already computed." << endl;
        Info <<
            "Check that the basis used for the parameterized BC are correct (RBF etc.)"
             << endl;
        Theta = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");

        M_Assert(Theta.cols() == Nbasis, "Reading wrong offline computations");
        for (label baseI = 0; baseI < Nbasis; baseI++)
        {
            Ttime.resize(0);
            ITHACAstream::read_fields(Ttime, "T_basis" + std::to_string(baseI + 1),
                                      basisFolderOffline);
            T_basis.append(Ttime.clone());
        }
    }
    else
    {
        Info << "\nComputing offline" << endl;
        Theta.resize(thermocouplesNum, Nbasis);
        offlineFlag = 1;
        timeSampleI = 0;

        Info << "Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;

        /// Tbasis
        for (label baseI = 0; baseI < Theta.cols(); baseI++)
        {
            Info << "\n--------------------------------------\n" << endl;
            Info << "Base " << baseI + 1 << " of " << Theta.cols() << endl;
            Info << "\n--------------------------------------\n" << endl;
            restart();
            heatFluxWeights = Foam::zero();
            heatFluxWeights[baseI] =  1;
            updateHeatFlux(heatFluxWeights);
            solveDirect();

            volScalarField spaceBase = list2Field(heatFlux[1]);
            for(int timeI = 0; timeI < offlineTimestepsSize; timeI++)
            {
                ITHACAstream::exportSolution(spaceBase,
                                             std::to_string(timeSteps[timeI + 1]),
                                             basisFolderOffline,
                                             "heatFluxBase" + std::to_string(baseI + 1));
                ITHACAstream::exportSolution(Ttime[timeI], 
                        std::to_string(timeSteps[timeI + 1]), basisFolderOffline, 
                        "T_basis" + std::to_string(baseI + 1));
            }
            T_basis.append(Ttime.clone());
        }

        for (label baseI = 0; baseI < Nbasis; baseI++)
        {
            Info << "Matrices: base " << baseI << endl;
            Eigen::VectorXd Tbase_vec =
                fieldValueAtThermocouples(T_basis[baseI][offlineTimestepsSize - 1]);
            for(int i = 0; i < thermocouplesNum; i++)
            {
                Theta(i, baseI) = Tbase_vec(i);
            }
        }

        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
    }
    M_Assert(T_basis.size() == Nbasis, "Somethong went wrong in the offline phase");
    offlineFlag = 0;
    Info << "\nOffline ENDED" << endl;
}

void sequentialIHTP_constant::reconstructT(word _outputFolder)
{
    Info << "Reconstructing field T" << endl;
    M_Assert(offlineFlag == 0, "Call during online phase");
    Ttime.resize(0);
    Info << "\nExporting solution in the time domain (" <<
        timeSteps[NtimeStepsBetweenSamples * timeSampleI] << ", " <<
        timeSteps[NtimeStepsBetweenSamples + NtimeStepsBetweenSamples * timeSampleI] <<
        "]\n" << endl;
    restart();
    if(timeSampleI == 0)
    {
        ITHACAstream::exportSolution(T_ic_time[0], std::to_string(timeSteps[0]),
                                     _outputFolder,
                                     reconstructedTemperatureName);
        volScalarField heatFlux_field = list2Field(heatFlux[0]);
        ITHACAstream::exportSolution(heatFlux_field,
                                     std::to_string(timeSteps[0]),
                                     _outputFolder,
                                     reconstructedHeatFluxName);
    }
    M_Assert(T_ic_ready == 1,
            "The T_ic field has not be computed for this iteration");

    for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
    {
        volScalarField T(_T);
        ITHACAutilities::assignIF(T, homogeneousBC);

        forAll(T_basis, baseI)
        {
            T += heatFluxTimeBasis[timeI][baseI] * T_basis[baseI][timeI];
        }
        T += T_ic_time[timeI];

        label realTimeStep = timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
        ITHACAstream::exportSolution(T, std::to_string(timeSteps[realTimeStep]),
                                     _outputFolder,
                                     reconstructedTemperatureName);
        volScalarField heatFlux_field = list2Field(heatFlux[realTimeStep]);
        ITHACAstream::exportSolution(heatFlux_field,
                                     std::to_string(timeSteps[realTimeStep]),
                                     _outputFolder,
                                     reconstructedHeatFluxName);
        Ttime.append(T.clone());
    }
    T_ic_ready = 0;
    Info << "ReconstructT END" << endl;
}

scalar sequentialIHTP_constant::reconstructTatThermocouple(label _TCindex)
{
    M_Assert(offlineFlag == 0, "Call during online phase");
    M_Assert(T_ic_ready == 1,
            "The T_ic field has not be computed for this iteration");

    Eigen::VectorXd weights = Foam2Eigen::List2EigenMatrix(heatFluxWeights);
    Eigen::VectorXd oldWeights = Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
    scalar out;
    Info << "Recontruction of T at thermocouple " << _TCindex + 1 << endl;
    Eigen::VectorXd T_ic_atTC =
        fieldValueAtThermocouples(T_ic_time[NtimeStepsBetweenSamples - 1]);
    Eigen::VectorXd Tvec = Theta * weights + T_ic_atTC; 

    out = Tvec(_TCindex);
    Info << "reconstructTatThermocouple END" << endl;
    return out;
}

void sequentialIHTP_constant::computeHeatFluxWeights(volScalarField _initialField, 
        List<scalar> _initialWeights, word _outputFolder)
{
    timeSampleI = 0;
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

    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        heatFluxWeightsOld = heatFluxWeights;

        if(timeSampleI == 0)
        {
            heatFluxWeightsGrad_old = heatFluxWeights * 0.0;
        }
        else
        {
            heatFluxWeightsGrad_old = heatFluxWeightsGrad;
        }

        solveT_ic(_initialField);

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
        Info << "Weights = \n" << heatFluxWeights << endl;
        updateHeatFlux(heatFluxWeights);

        heatFluxWeightsGrad = (heatFluxWeights - heatFluxWeightsOld) /
            timeSamplesDeltaT;

        reconstructT(_outputFolder);

        label verbose = 0;
        parameterizedHeatFlux_postProcess(linSys, weigths, _initialField, _outputFolder, 
                verbose);
        timeSampleI++;
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms =
            std::chrono::duration<double, std::milli>(t_end-t_start).count();
        Info << "CPU time = " << elapsed_time_ms << " milliseconds" << endl << endl;
    }
    Info << "End" << endl;
    Info << endl;
}

void sequentialIHTP_constant::solveT_ic(volScalarField _initialField)
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
