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
            heatFluxSpaceBasis[funcI][faceI] = 1e6 * Foam::sqrt(1 + (shapeParameter *
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
            forAll(tempBasis[baseI], cellI)
            {
                tempBasis[baseI][cellI] *= 1e6;
            }
            volScalarField base = list2Field(tempBasis[baseI], 0.0);
            ITHACAstream::exportSolution(base,
                                         std::to_string(1),
                                         debugFolder,
                                         "PODbase" + std::to_string(baseI + 1));
        }
        heatFluxSpaceBasis = tempBasis;
    }
}

void sequentialIHTP::set_gParametrized(word spaceBaseFuncType,
        scalar shapeParameter_space, label Npod)
{
    volScalarField& T = _T();
    setSpaceBasis(spaceBaseFuncType, shapeParameter_space, Npod);

    NsamplesWindow = 1;
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

List<List<scalar>> sequentialIHTP::interpolateWeights(List<scalar> Wold, List<scalar> Wnew)
{
    M_Assert(Wold.size() == Wnew.size(), "Input weights vectors must have the same size");

    double t0 = 0;
    double t1 = NtimeStepsBetweenSamples * deltaTime;
    List<List<scalar>> Wout;
    Wout.resize(Wold.size());
    forAll (Wold, wI)
    {
        Wout[wI].resize(NtimeStepsBetweenSamples + 1);
        for(int timeI = 0; timeI < NtimeStepsBetweenSamples + 1; timeI++)
        {
            double time = (timeI + 1) * deltaTime;
            double a = Wold[wI] - (Wnew[wI] - Wold[wI]) / (t1 - t0) * t0;
            double b = (Wnew[wI] - Wold[wI]) / (t1 - t0);
            Wout[wI][timeI] = a + b * time;
        }
    }
    return Wout;
}

void sequentialIHTP::updateTimeHeatFlux(List<scalar> weights)
{
    if(offlineFlag)
    {
        forAll (weights, weightI)
        {
            for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
            {
                heatFluxTimeBasis[timeI][weightI] = weights[weightI]; 
            }
        }
    }
    else
    {
        scalar oldSamplingTime; 
        if(timeSampleI == 0)
        {
            oldSamplingTime = timeSteps[0];
        }
        else
        {
            oldSamplingTime = samplingTime[timeSampleI - 1];
        }
        Info << "debug: oldSamplingTime = " << oldSamplingTime << endl;
        forAll (weights, weightI)
        {
            for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
            {
                if(linearBasis)
                {
                    label realTimeStep;
                    if(timeSampleI == 0)
                    {
                        realTimeStep =  timeI + 1;
                    }
                    else
                    {
                        realTimeStep = samplingSteps[timeSampleI - 1] + timeI + 1;
                    }
                    scalar realTime = timeSteps[realTimeStep];
                    //Info << "debug: realTime = " << realTime << endl;
                    heatFluxTimeBasis[timeI][weightI] = heatFluxWeightsOld[weightI] + 
                        (realTime - oldSamplingTime) * (weights[weightI] - 
                                heatFluxWeightsOld[weightI]) / timeSamplesDeltaT; 
                }
                else // Piecewise constant heatflux in time
                {
                    heatFluxTimeBasis[timeI][weightI] = weights[weightI]; 
                }
            }
        }
    }
}

void sequentialIHTP::updateHeatFlux(List<scalar> weights)
{
    M_Assert(weights.size() == Nbasis,
             "weigths size different from basis functions size");
    volScalarField& T = _T();

    label firstTimeI = timeSampleI * NtimeStepsBetweenSamples;
    int lastTimeStep = firstTimeI + NtimeStepsBetweenSamples;

    updateTimeHeatFlux(weights);
    for(int timeI = firstTimeI; timeI < lastTimeStep; timeI++)
    {    
        //Info << "debug: timeI = " << timeI << endl;
        //Info << "debug: firstTimeI = " << firstTimeI << endl;
        //Info << "debug: timeI - firstTimeI = " << timeI - firstTimeI << endl;
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
        //forAll (weights, weightI)
        //{
        //    Info << "debug: heatFluxTimeBasis[" << timeI - firstTimeI << "][" << weightI << "] = " << heatFluxTimeBasis[timeI - firstTimeI][weightI] << endl;
        //}
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

void sequentialIHTP::parameterizedBCoffline(bool force)
{
    fvMesh& mesh = _mesh();
    Tbasis.resize(Nbasis);
    Ttau.resize(0);
    T0field.resize(0);
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
        Theta_tau = ITHACAstream::readMatrix(folderOffline + "Theta_tau_mat.txt");
        PtrList<volScalarField> Ttemp;
        ITHACAstream::read_fields(Ttemp, "T", basisFolderOffline);
        M_Assert(Ttemp.size() == Nbasis, 
                "The Tbasis read from file has not the right size");

        M_Assert(Theta.cols() == Nbasis, "Reading wrong offline computations");
        for (label baseI = 0; baseI < Nbasis; baseI++)
        {
            Tbasis.set(baseI, Ttemp[baseI].clone());
            Ttime.resize(0);
            ITHACAstream::read_fields(Ttime, "Ttau" + std::to_string(baseI + 1),
                                      folderOffline);
            Ttau.append(Ttime.clone());
        }
    }
    else
    {
        Info << "\nComputing offline" << endl;
        Theta.resize(thermocouplesNum, Nbasis);
        Theta_tau = Theta;
	offlineFlag = 1;
	timeSampleI = 0;

        Info << "Theta size = " << Theta.rows() << ", " << Theta.cols() << endl;

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
            
            volScalarField& T = _T();
            Tbasis.set(baseI, T.clone());
            volScalarField gParametrizedField = list2Field(heatFlux[1]);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(baseI + 1),
                                         basisFolderOffline,
                                         "g");
            ITHACAstream::exportSolution(T, std::to_string(baseI + 1),
                                         basisFolderOffline,
                                         "T");
        }

        for (label baseI = 0; baseI < Theta.cols(); baseI++)
        {
            Info << "\n--------------------------------------\n" << endl;
            Info << "Base " << baseI + 1 << " of " << Theta.cols() << endl;
            Info << "\n--------------------------------------\n" << endl;
            restart();
            Ttime.resize(0);
            heatFluxWeights = Foam::zero();
            
            // Compute Ttau
            solveTtau(baseI);
            M_Assert(Ttime.size() == offlineTimestepsSize, "Wrong restert time for Ttau");
            M_Assert(Ttau_lastTime.size() == 1, "Wrong restert time for Ttau");
            for(int timeI = 0; timeI < offlineTimestepsSize; timeI++)
            {
                volScalarField& T = Ttime[timeI];
                ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI + 1]),
                                             folderOffline,
                                             "Ttau" + std::to_string(baseI + 1));
            }
            Ttau.append(Ttime.clone());
            Tcomp = fieldValueAtThermocouples(Ttime[offlineTimestepsSize - 1]);
            for(int i = 0; i < Tcomp.size(); i++)
            {
                Theta_tau(i, baseI) = 1.0 / timeSamplesDeltaT * Tcomp(i);
            }
            Tcomp = fieldValueAtThermocouples(Tbasis[baseI]);
            for(int i = 0; i < Tcomp.size(); i++)
            {
                Theta(i, baseI) = Tcomp(i) + Theta_tau(i, baseI);
            }
        }

        ITHACAstream::exportMatrix(Theta, "Theta", "eigen", folderOffline);
        ITHACAstream::exportMatrix(Theta_tau, "Theta_tau", "eigen", folderOffline);
        Eigen::MatrixXd ThetaTheta_tau = Theta + Theta_tau;
        ITHACAstream::exportMatrix(ThetaTheta_tau, "ThetaTheta_tau", "eigen", 
                folderOffline);
    }

    Eigen::MatrixXd A = Theta;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singularValues = svd.singularValues();
    double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
    Info << "Theta Condition number = " << conditionNumber << endl;
    ITHACAstream::exportMatrix(singularValues, "ThetaSingularValues", "eigen",
                               folderOffline);
    offlineFlag = 0;
    Info << "\nOffline ENDED" << endl;
}

void sequentialIHTP::reconstructTF()
{
    TF.resize(0);
    for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
    {
        volScalarField T(_T);
        ITHACAutilities::assignIF(T, homogeneousBC);
        forAll(heatFluxWeights, baseI)
        {
            scalar coeff = heatFluxWeights[baseI] - heatFluxWeightsOld[baseI];
            coeff = coeff / timeSamplesDeltaT; 
            T += coeff * Ttau[baseI][timeI];
        }
        TF.append(T.clone());
    }
}

void sequentialIHTP::solveTF(word _outputFolder, volScalarField _initialField)
{
    Info << "Solving TF" << endl;
    volScalarField TF_old = _initialField;

    // Set the first element because it will export it in reconstructT
    T0_time.append(_initialField.clone());
    
    if(timeSampleI > 0)
    {
        TF_old = TF[TF.size() - 1]; 
    }
    TF.resize(0);
    restartOffline();
    M_Assert(diffusivity > 1e-36, "Set the diffusivity value");
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    dimensionedScalar dt("dt", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0);

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
    if(timeSampleI == 0)
    {
        forAll(Tbasis, baseI)
        {
            _initialField -= heatFluxWeightsOld[baseI] * Tbasis[baseI];
        }
        ITHACAutilities::assignIF(T, _initialField);
    }
    else
    {
        ITHACAutilities::assignIF(T, TF_old);
    }

    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();
    fv::options& fvOptions(_fvOptions());
    label timeI = 0;

    volScalarField source = dt * T;
    ITHACAutilities::assignIF(source, homogeneousBC);
    forAll(Tbasis, baseI)
    {
        scalar coeff = heatFluxWeights[baseI] - heatFluxWeightsOld[baseI];
        coeff = coeff / timeSamplesDeltaT; 
        source += coeff * dt * Tbasis[baseI];
    }

    while (runTime.loop())
    {
        Info << "\nTime = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT * diffusivity, T)
                ==
                - source 
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        label realTimeStep = timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
        ITHACAstream::exportSolution(T, std::to_string(timeSteps[realTimeStep]),
                                     _outputFolder,
                                     "TF");
        ITHACAstream::exportSolution(source, std::to_string(timeSteps[realTimeStep]),
                                     _outputFolder,
                                     "TFsource");
        ITHACAstream::exportSolution(_initialField, std::to_string(timeSteps[realTimeStep]),
                                     _outputFolder,
                                     "initialField");
        TF.append(T.clone());

        runTime.printExecutionTime(Info);
        runTime.write();
        timeI++;
    }
    TF_ready = 1;
    Info << "TF computation ENDED" << endl << endl;
}

void sequentialIHTP::reconstrucT(word outputFolder)
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
        ITHACAstream::exportSolution(T0_time[0], std::to_string(timeSteps[0]),
                                     outputFolder,
                                     "Treconstructed");
        volScalarField gParametrizedField = list2Field(heatFlux[0]);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(timeSteps[0]),
                                     outputFolder,
                                     "gReconstructed");
    }
    M_Assert(TF_ready == 1, 
            "The TF field has not be computed or reconstructed for this iteration");

    for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
    {
        volScalarField T(_T);
        ITHACAutilities::assignIF(T, homogeneousBC);
        if(linearBasis == 1)
        {
            Info << "LINEAR recontruction of T" << endl;
            forAll(Tbasis, baseI)
            {
                T += heatFluxTimeBasis[timeI][baseI] * Tbasis[baseI]; 
            }
            T += TF[timeI];
        }
        else
        {
            Info << "CONSTANT basis not yet implemented, exiting" << endl;
            exit(78);
            
        }

        label realTimeStep = timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
        ITHACAstream::exportSolution(T, std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "Treconstructed");
        volScalarField gParametrizedField = list2Field(heatFlux[realTimeStep]);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(timeSteps[realTimeStep]),
                                     outputFolder,
                                     "gReconstructed");
        Ttime.append(T.clone());
    }
    TF_ready = 0;
}

scalar sequentialIHTP::reconstrucT(vector _point, label _timeI)
{
    label realTimeStep = _timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
    Info << "Reconstructing field T at point " << _point << " and time " << 
        timeSteps[realTimeStep] << endl;
    M_Assert(offlineFlag == 0, "Call during online phase");

    scalar out = 0;
    scalar TF_point = fieldValueAtPoint(TF[_timeI], _point);
    if(linearBasis == 1)
    {
        Info << "LINEAR recontruction of T" << endl;
        forAll(Tbasis, baseI)
        {
            scalar Tbasis_point = fieldValueAtPoint(Tbasis[baseI], _point);
            out += heatFluxTimeBasis[_timeI][baseI] * Tbasis_point;
        }
        out += TF_point;
    }
    else
    {
        Info << "CONSTANT basis not yet implemented, exiting" << endl;
        exit(78);
    }
    Info << "DONE" << endl;
    return out;
}

scalar sequentialIHTP::reconstrucT(vector _point)
{
    label _timeI = NtimeStepsBetweenSamples - 1;
    label realTimeStep = _timeI + NtimeStepsBetweenSamples * timeSampleI + 1;
    Info << "Reconstructing field T at point " << _point << " and time " << 
        timeSteps[realTimeStep] << endl;
    M_Assert(offlineFlag == 0, "Call during online phase");

    scalar out = 0;
    Eigen::VectorXd ThetaPrime(Nbasis);
    Eigen::VectorXd ThetaTauPrime(Nbasis);
    Eigen::VectorXd ThetaTildePrime(Nbasis);
    Eigen::VectorXd weights = Foam2Eigen::List2EigenMatrix(heatFluxWeights);
    Eigen::VectorXd oldWeights = Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
    forAll(Tbasis, baseI)
    {
        ThetaPrime(baseI) = fieldValueAtPoint(Tbasis[baseI], _point);
        ThetaTauPrime(baseI) = fieldValueAtPoint(Ttau[baseI][_timeI], _point);
        ThetaTildePrime(baseI) = ThetaPrime(baseI) + 
            (1.0 / timeSamplesDeltaT) * ThetaTauPrime(baseI);
    }
    if(linearBasis == 1)
    {
        Info << "LINEAR recontruction of T" << endl;
        scalar temp = (1.0 / timeSamplesDeltaT);
        scalar a = ThetaTildePrime.dot(weights);
        scalar b = ThetaTauPrime.dot(oldWeights);
        out = a - temp * b;
    }
    else
    {
        Info << "CONSTANT basis not yet implemented, exiting" << endl;
        exit(78);
    }
    Info << "DONE" << endl;
    return out;
}

scalar sequentialIHTP::reconstrucTatThermocouple(label _TCindex)
{
    M_Assert(offlineFlag == 0, "Call during online phase");

    Eigen::VectorXd weights = Foam2Eigen::List2EigenMatrix(heatFluxWeights);
    Eigen::VectorXd oldWeights = Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
    scalar out;
    if(linearBasis == 1)
    {
        Info << "LINEAR recontruction of T" << endl;
        Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
        Eigen::MatrixXd ThetaTilde = Theta + Temp;
        Eigen::VectorXd Tvec = ThetaTilde * weights + T0_vector - 
        Temp * oldWeights;
        out = Tvec(_TCindex);
    }
    else
    {
        Info << "CONSTANT basis not yet implemented, exiting" << endl;
        exit(78);
    }
    Info << "DONE" << endl;
    return out;
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

Eigen::VectorXd sequentialIHTP::reconstrucT(Eigen::VectorXi cells)
{
    Eigen::VectorXd Tout(cells.size());
    int timeI = NtimeStepsBetweenSamples - 1;
    for(int cellI = 0; cellI < cells.size(); cellI++)
    {
        forAll(Tbasis, baseI)
        {
            Tout(cellI) += heatFluxWeights[baseI] * 
                Tbasis[baseI].internalField()[cellI];
        }
        Tout(cellI) += TF[timeI].internalField()[cellI];
    }

    return Tout;
}

void sequentialIHTP::parameterizedBC(word outputFolder, volScalarField _initialField,
        List<scalar> _initialWeights)
{
    Info << endl << "Using quasilinearity of direct problem ::" << endl;
    Info << "Using " << linSys_solver << " to solve the linear system" << endl;
    if(linearBasis == 0)
    {
        Info << "\nCONSTANT time basis\n" << endl;
        computeConstantHeatWeights(outputFolder, _initialField);
    }
    else
    {
        Info << "\nLINEAR time basis\n" << endl;
        computeLinearHeatWeights(outputFolder, _initialField, _initialWeights);
    }
}

Eigen::VectorXd sequentialIHTP::solveLinSys(List<Eigen::MatrixXd> linSys)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
    ITHACAregularization::Picard(svd.matrixU(), svd.singularValues(),
            linSys[1]);
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
        linSys = ITHACAregularization::precondition(linSys, "SVD");

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

void sequentialIHTP::computeConstantHeatWeights(word outputFolder, 
        volScalarField _initialField)
{
    M_Assert(linearBasis == 0, "Wrong call of function computeConstantHeatWeights");
    timeSampleI = 0;

    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    linSys[0] = Theta.transpose() * Theta;
    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        auto t_start = std::chrono::high_resolution_clock::now();


        if(timeSampleI > 0)
        {
             ITHACAutilities::assignIF(_initialField, Ttime[NtimeStepsBetweenSamples -1]);
        }
	solveT0(_initialField);

	TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, thermocouplesNum * 
                NsamplesWindow);
        linSys[1] = Theta.transpose() * (TmeasShort - T0_vector);

        ITHACAstream::exportMatrix(linSys[0], "linSys0." + std::to_string(timeSampleI), 
                "eigen", outputFolder);
        ITHACAstream::exportMatrix(linSys[1], "linSys1." + std::to_string(timeSampleI), 
                "eigen", outputFolder);

        Eigen::VectorXd weigths = solveLinSys(linSys);

        ITHACAstream::exportMatrix(weigths, "weigths" + std::to_string(timeSampleI), 
                "eigen", outputFolder);

        heatFluxWeights.resize(weigths.size());
        forAll(heatFluxWeights, weightI)
        {
            heatFluxWeights[weightI] = weigths(weightI);
        }
        Info << "Weights = \n" << heatFluxWeights << endl;
        updateHeatFlux(heatFluxWeights);
        label verbose = 0;
        parameterizedBC_postProcess(linSys, weigths, _initialField, outputFolder, verbose);
	timeSampleI++;
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms = 
            std::chrono::duration<double, std::milli>(t_end-t_start).count();
        Info << "CPU time = " << elapsed_time_ms << " milliseconds" << endl << endl;
    }
    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", outputFolder);
    Info << "End" << endl;
    Info << endl;
}

void sequentialIHTP::computeLinearHeatWeights(word outputFolder, 
        volScalarField _initialField, List<scalar> _initialWeights)
{
    M_Assert(linearBasis == 1, "Wrong call of function computeLinearHeatWeights");
    M_Assert(_initialWeights.size() == Nbasis, "Wrong input weights");
    timeSampleI = 0;

    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    heatFluxWeightsOld = _initialWeights;

    while(timeSampleI < timeSamplesNum)
    {
        Info << "\nTime sample " << timeSampleI + 1 << endl;
        auto t_start = std::chrono::high_resolution_clock::now();

        if(timeSampleI > 0)
        {
             ITHACAutilities::assignIF(_initialField, Ttime[NtimeStepsBetweenSamples -1]);
             heatFluxWeightsOld = heatFluxWeights;
        }
	solveT0(_initialField);

	TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, thermocouplesNum * 
                NsamplesWindow);
        Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
        Eigen::MatrixXd M = Theta + Temp;
        Eigen::VectorXd heatFluxWeightsOld_Eig = 
            Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
        linSys[0] = M.transpose() * M;
        linSys[1] = M.transpose() * (TmeasShort + 
                Temp * heatFluxWeightsOld_Eig - T0_vector);
        ITHACAstream::exportMatrix(linSys[0], "linSys0." + std::to_string(timeSampleI), 
                "eigen", outputFolder);
        ITHACAstream::exportMatrix(linSys[1], "linSys1." + std::to_string(timeSampleI), 
                "eigen", outputFolder);
        ITHACAstream::exportMatrix(heatFluxWeightsOld_Eig, "heatFluxWeightsOld" + 
                std::to_string(timeSampleI), "eigen", outputFolder);

        Eigen::VectorXd weigths = solveLinSys(linSys);

        ITHACAstream::exportMatrix(weigths, "weigths" + std::to_string(timeSampleI), 
                "eigen", outputFolder);

        heatFluxWeights.resize(weigths.size());
        forAll(heatFluxWeights, weightI)
        {
            heatFluxWeights[weightI] = weigths(weightI);
        }
        Info << "Weights = \n" << heatFluxWeights << endl;
        updateHeatFlux(heatFluxWeights);
        label verbose = 0;
        parameterizedBC_postProcess(linSys, weigths, _initialField, outputFolder, verbose);
	timeSampleI++;
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_time_ms = 
            std::chrono::duration<double, std::milli>(t_end-t_start).count();
        Info << "CPU time = " << elapsed_time_ms << " milliseconds" << endl << endl;
    }
    ITHACAstream::exportMatrix(Jlist, "costFunction", "eigen", outputFolder);
    Info << "End" << endl;
    Info << endl;
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

void sequentialIHTP::solveT0(volScalarField _initialField)
{
    Info << "\nSolving FULL T0 problem" << endl;
    restartOffline();
    fvMesh& mesh = _mesh();
    simpleControl& simple = _simple();
    fv::options& fvOptions(_fvOptions());
    volScalarField T0(_T);
    Foam::Time& runTime = _runTime();
    set_valueFraction();
    List<scalar> RobinBC = Tf * 0.0;
    word outputFolder = "./ITHACAoutput/debugT0/";

    ITHACAutilities::assignIF(T0, _initialField);

    T0field.append(T0.clone());
    T0_time.resize(0);
    label timeI = 0;
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

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        timeI++;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T0) - fvm::laplacian(DT * diffusivity, T0)
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T0);
        }

        T0_time.append(T0.clone());
        T0field.append(T0.clone());
        //ITHACAstream::exportSolution(T0, std::to_string(
        //            timeSteps[samplingSteps[timeSampleI] - NtimeStepsBetweenSamples + 
        //            timeI]), outputFolder, "T0");
        runTime.printExecutionTime(Info);
        runTime.write();
    }

    T0_vector = fieldValueAtThermocouples(T0_time);
    Info << "SolveT0 ENDED\n" << endl;
}

void sequentialIHTP::getT0modes()
{
    word outputFolder = "./ITHACAoutput/modes/";
    Info << "Computing " << NmodesT0 << " T0 modes" << endl;

    ITHACAPOD::getModes(T0field, T0modes, "T0",
                        0, 0, 0,
                        NmodesT0);
    PtrList<volScalarField> modes = T0modes.toPtrList();
    forAll(modes, modeI)
    {
        ITHACAstream::exportSolution(modes[modeI], std::to_string(modeI + 1),
                outputFolder, "T0modes");
    }
    Info << "T0 modes COMPUTED\n" << endl;
}

void sequentialIHTP::projectT0()
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

    T0modes.toEigen();
    T0implicitMatrix_red = T0modes.EigenModes[0].transpose() 
        * T0implicitMatrix * T0modes.EigenModes[0];
    T0explicitMatrix_red = T0modes.EigenModes[0].transpose() 
        * b.asDiagonal() * T0modes.EigenModes[0];
    Info << "Projection matrices COMPUTED\n" << endl;
}

void sequentialIHTP::projectDirectOntoT0()
{
    /// Creation of the matrices to project direct solution at the last timestep
    /// onto the T0 reduced space
    Info << "Computing direct problem projection matrices" << endl;
    int internalFieldSize = Tbasis[0].internalField().size();
    Eigen::MatrixXd Tbasis_Eigen(internalFieldSize, Nbasis);

    PtrList<volScalarField> Tbasis_lastTime;
    M_Assert(Tbasis[0].size() == NtimeStepsBetweenSamples, 
            "The basis for the direct problem have wrong dimention in time");

    forAll(Tbasis, baseI)
    {
        volScalarField temp = Tbasis[baseI]; 
        Tbasis_lastTime.append(temp.clone());
    }

    Tbasis_projectionMat = T0modes.project(Tbasis_lastTime); 
}

void sequentialIHTP::pointProjectionOffline()
{
    int Ncells = magicPoints.size();
    M_Assert(Ncells > 0, "Set the number of magic points");

    int lastTimestepID = NtimeStepsBetweenSamples - 1;

    pointsReconstructMatrix.resize(Ncells, NmodesT0);
    for(int cellI = 0; cellI < Ncells; cellI++)
    {
        pointsReconstructMatrix.row(cellI) = 
            T0modes.EigenModes[0].row(magicPoints[cellI]);
    }
    pointTbasis_reconstructionMat = pointsReconstructMatrix * Tbasis_projectionMat;
    pointTad_reconstructed = pointsReconstructMatrix * Tad_projected;
}

void sequentialIHTP::projectionErrorOffline()
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
        volScalarField base = Tbasis[baseI];
        volScalarField baseProj(base);
        T0modes.projectSnapshot(base, baseProj, NmodesT0, "L2");
        volScalarField temp = base - baseProj;
        projectionErrorTbasis.append(temp.clone());
        ITHACAstream::exportSolution(projectionErrorTbasis[baseI], std::to_string(baseI + 1),
                outputFolder, "projectionErrorTbasis");
    }
    projectionErrorTad.resize(0);
    volScalarField TadProj(Tbasis[0]);
}

void sequentialIHTP::T0offline(int NmagicPoints)
{
    getT0modes();
    findMagicPoints(NmagicPoints);
    projectT0();
    projectDirectOntoT0();
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

    dimensionedScalar Cond("Cond", dimensionSet(1, 1, -3, -1, 0, 0, 0), thermalCond);

    while (simple.loop())
    {
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::laplacian(Cond, T)
            );
            TEqn.solve();
        }
    }
    Info << "Direct computation ENDED" << endl;
    
}

void sequentialIHTP::solveTtau(label _baseI)
{
    Info << "Solving Ttau" << endl;
    M_Assert(offlineFlag, "solveTtau should be called only during offline phase");
    
    //restartOffline();
    restartOfflineTau();
    M_Assert(diffusivity > 1e-36, "Set the diffusivity value");
    volScalarField& T = _T();
    fvMesh& mesh = _mesh();
    dimensionedScalar dt("dt", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0);
    List<scalar> RobinBC = Tf * 0.0;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
        {
            ITHACAutilities::assignMixedBC(T, patchI, RobinBC, refGrad,
                                           valueFraction);
        }
        else
        {
            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
        }
    }
    ITHACAutilities::assignIF(T, homogeneousBC);

    simpleControl& simple = _simple();
    Foam::Time& runTime = _runTime();
    fv::options& fvOptions(_fvOptions());
    label timeI = 0;
    Ttime.resize(0);
    Ttau_lastTime.resize(0);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        volScalarField source = dt * T;
        ITHACAutilities::assignIF(source, homogeneousBC);
        if(timeI > 0)
        {
            source = dt * Tbasis[_baseI];
        }

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT * diffusivity, T)
                ==
                - source 
            );
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        if(timeI <= offlineTimestepsSize - 1)
        {
            Ttime.append(T.clone());
        }
        else if(timeI == offlineTimestepsSize)
        {
            Ttau_lastTime.append(T.clone());
        }
        else
        {
            Info << "Wrong restart of Ttau, EXITING" << endl;
            exit(10);
        }

        runTime.printExecutionTime(Info);
        runTime.write();
        timeI++;
    }
    Info << "Ttau computation ENDED" << endl;
}

void sequentialIHTP::readThermocouples()
{
    Info << "Defining positions of thermocouples" << endl;

    if (!thermocouplesRead)
    {
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
            fieldInt.segment(sampleTimeI * thermocouplesNum, thermocouplesNum) =
                fieldValueAtThermocouples(fieldList, samplingSteps[sampleTimeI]);
        }
    }
    else if ( fieldList.size() == NtimeStepsBetweenSamples )
    {
        Info << "\nField size = " << fieldList.size() << ".\nSampling ONLY the last timestep\n\n" << endl;
        fieldInt =
            fieldValueAtThermocouples(fieldList, NtimeStepsBetweenSamples - 1);
    }
    else
    {
        Info << "The input fieldList of sequentialIHTP::fieldValueAtThermocouples can have size Ntimes + 1 (=" << Ntimes + 1 << ") or\n";
        Info << " NtimeStepsBetweenSamples  (=" <<  NtimeStepsBetweenSamples << ") but has size " << fieldList.size() << endl;
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

void sequentialIHTP::restartOfflineTau()
{
    Info << "Setting endTime to offlineEndTime + 1" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(offlineEndTime + deltaTime);
    Info << "Ready for new offline computation" << endl;
}

void sequentialIHTP::restartT0()
{
    Info << "Setting endTime to offlineEndTime" << endl;
    restart();
    Time& runTime = _runTime();
    instantList Times = runTime.times();
    runTime.setTime(0.0, 0);
    runTime.setEndTime(timeSamplesDeltaT);
    Info << "Ready for new T0 computation" << endl;
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

void sequentialIHTP::parameterizedBC_postProcess(
    List<Eigen::MatrixXd> linSys, Eigen::VectorXd weigths, volScalarField _initialField, 
    word outputFolder, label verbose)
{
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
        std::cout << "T0_vector = " << std::endl;
        std::cout << T0_vector << std::endl;
        std::cout << "Tmeas = " << std::endl;
        std::cout << Tmeas << std::endl;
    }

    reconstrucT(outputFolder);
    Tcomp = fieldValueAtThermocouples(Ttime);
    std::cout << "Tcomp = \n" << Tcomp.transpose() << std::endl;
    std::cout << "TmeasShort = \n" << TmeasShort.transpose() << std::endl;
    J = 0.5 * Foam::sqrt((Tcomp - TmeasShort).dot(Tcomp - TmeasShort));
    Info << "J = " << J << endl;
    Jlist.conservativeResize(Jlist.size() + 1);
    Jlist(Jlist.size() - 1) = J;
}

void sequentialIHTP::findMagicPoints(int NmagicPoints)
{
    M_Assert(NmagicPoints > 0, "Set number of magic points");
    magicPoints.clear();
    M_Assert(NmagicPoints <= NmodesT0, 
            "Number of magic points bigger than number of modes");
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd r;
    Eigen::VectorXd rho(1);
    Eigen::MatrixXd MatrixModes = T0modes.toEigen()[0];
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
