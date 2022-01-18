#ifndef seqTest_H
#define seqTest_H

#include "ReducedSequentialIHTP_incrementalPOD.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"

class seqTest: public reducedSequentialIHTP_incrementalPOD
{
    public:
        explicit seqTest(int argc, char* argv[])
            :
            reducedSequentialIHTP_incrementalPOD(argc, argv),
            initialField(_T())
        {
            fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
            folderOffline = "./ITHACAoutput/offlineParamBC_unsteady/";
        }
        word directOutputFolder = "./ITHACAoutput/directProblem/";
        scalar a = 0;
        scalar b = 0;
        scalar c = 0;
        scalar d = 0;
        PtrList<volScalarField> Ttrue;
        PtrList<volScalarField> recT;
        PtrList<volScalarField> solT;
        volScalarField initialField;

        Foam::vector probeCenter = Foam::vector(1, 0.0, 0.6);
        Foam::vector probeSide = Foam::vector(0.8, 0.0, 1.1);
        Foam::vector probe2 = Foam::vector(1.0, 0.02, 0.6);

        Eigen::MatrixXd probeCenter_trueHeatFlux;
        Eigen::MatrixXd probeCenter_gRec;
        Eigen::MatrixXd probeCenter_Ttrue;
        Eigen::MatrixXd probeCenter_Trec;

        Eigen::MatrixXd probeSide_trueHeatFlux;
        Eigen::MatrixXd probeSide_gRec;

        List<scalar> heatFlux_steady;
        scalar maxFrequency; // Hz

        Eigen::MatrixXi List2EigenMatrix ( List<int> list )
        {
            Eigen::MatrixXi matrix(list.size(), 1);

            for (label i = 0; i < matrix.rows(); i++)
            {
                matrix(i, 0) = list[i];
            }

            return matrix;
        }

        Eigen::MatrixXd List2EigenMatrix ( List<scalar> list )
        {
            Eigen::MatrixXd matrix(list.size(), 1);

            for (label i = 0; i < matrix.rows(); i++)
            {
                matrix(i, 0) = list[i];
            }

            return matrix;
        }

        void assignTrueIF()
        {
            volScalarField& T(_T());

            for (label i = 0; i < T.internalField().size(); i++)
            {
                T.ref()[i] =  initialField.internalField()[i];
            }
        }

        /// Sets steady component of heat flux, g
        void set_steadyHeatFlux()
        {
            fvMesh& mesh(_mesh());
            volScalarField& T(_T());
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
            heatFlux_steady.resize(hotSideSize);
            forAll(heatFlux_steady, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                scalar faceZ =
                    mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                heatFlux_steady[faceI] = 
                    - thermalCond * (a * faceX * faceX + b * faceZ + c);
            }
        }

        void set_trueHeatFlux(scalar timeGrad)
        {
            fvMesh& mesh(_mesh());
            volScalarField& T(_T());
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
            trueHeatFlux.resize(timeSteps.size());
            set_steadyHeatFlux();
            forAll(timeSteps, timeI)
            {
                // Setting trueHeatFlux
                trueHeatFlux[timeI].resize(hotSideSize);
                forAll(trueHeatFlux[timeI], faceI)
                {
                    scalar faceX =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                    scalar faceZ =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();

                    //scalar g1 = - thermalCond * (b * faceZ * faceZ  + c);
                    //trueHeatFlux[timeI][faceI] = timeGrad * timeSteps[timeI] * g1 + g1; 

                    scalar g1 = - thermalCond * (b * faceZ * faceZ  + c);
                    scalar g2 = - thermalCond * c / (1 + (faceX - 1) * (faceX - 1) + 
                            faceZ * faceZ);
                    trueHeatFlux[timeI][faceI] = g1 + 
                        g1 / 2 * (std::sin ( 2 * M_PI * maxFrequency * 
                                    timeSteps[timeI] * timeSteps[timeI] / endTime)) +
                        g2 * std::exp(- 0.1 * timeSteps[timeI]);
                }
            }
        }

        void set_Tf(label timeI)
        {
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());
            // Setting BC at the cold side
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            label coldSideSize = T.boundaryField()[coldSide_ind].size();
            Tf.resize(coldSideSize);
            forAll(Tf, faceI)
            {
                Tf[faceI] = 300;
            }
            set_valueFraction();
        }

        void assignDirectBC(label timeI)
        {
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());
            set_Tf(timeI);
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
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

        void solveTrue()
        {
            Info << "Performing true computation" << endl;
            restart();
            volScalarField& T(_T());
            Foam::Time& runTime = _runTime();
            Info << "Thermal diffusivity = " << diffusivity << " m2/s" << endl;
            simpleControl& simple = _simple();
            fv::options& fvOptions(_fvOptions());
            label timeI = 0;
            assignTrueIF();
            heatFlux = trueHeatFlux;
            assignDirectBC(timeI);
            ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI]),
                                         directOutputFolder,
                                         "Ttrue");
            volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
            ITHACAstream::exportSolution(trueHeatFluxField,
                                         std::to_string(timeSteps[timeI]), 
                                         directOutputFolder,
                                         "trueHeatFlux");
            Ttrue.resize(0);
            Ttrue.append(T.clone());

            while (runTime.loop())
            {
                Info << "Time = " << runTime.timeName() << nl << endl;
                timeI++;
                assignDirectBC(timeI);

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

                ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI]),
                                             directOutputFolder,
                                             "Ttrue");
                trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                ITHACAstream::exportSolution(trueHeatFluxField,
                                             std::to_string(timeSteps[timeI]), 
                                             directOutputFolder,
                                             "trueHeatFlux");
                Ttrue.append(T.clone());
                Info << endl;
                runTime.printExecutionTime(Info);
                runTime.write();
            }
            

            Tmeas = fieldValueAtThermocouples(Ttrue);
            Info << "End\n" << endl;
        }

        void addNoise(double noiseStdDev_)
        {
            double variance = Tmeas.mean() * noiseStdDev_ * noiseStdDev_;
            std::cout << "debug : noiseStdDev_ = " << noiseStdDev_  << std::endl;
            std::cout << "debug : Tmeas.mean() = " << Tmeas.mean()  << std::endl;
            std::cout << "debug : variance = " << variance  << std::endl;
            int Ntc = thermocouplesPos.size();

            auto noiseDensity = std::make_shared<muq::Modeling::Gaussian>
                                (Eigen::VectorXd::Zero(Ntc), 
                                 Eigen::MatrixXd::Identity(Ntc, Ntc) * variance);
            forAll(samplingSteps, tI)
            {
                Tmeas.segment(tI * Ntc, Ntc) += noiseDensity->Sample();
            }
        }

        void inverseProblemPostProcess(word outputFolder)
        {
            fvMesh& mesh(_mesh());
            Eigen::VectorXd heatFluxL2norm(timeSteps.size());
            Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());
            Eigen::VectorXd trueHeatFluxL2norm(timeSteps.size());
            Eigen::VectorXd trueHeatFluxLinfNorm(timeSteps.size());
            Eigen::VectorXd avgRelErrL2norm(timeSteps.size());
            Eigen::VectorXd avgRelErrLinfNorm(timeSteps.size());
            Info << "debug : g.size() = " << heatFlux.size() << endl;
            Info << "debug : trueHeatFlux.size() = " << trueHeatFlux.size() << endl;
            Info << "debug : Ttrue.size() = " << Ttrue.size() << endl;
            Info << "debug : timeSteps.size() = " << timeSteps.size() << endl;

            label cellprobeCenter = mesh.findCell(probeCenter);
            Eigen::VectorXd probeCenter_trueHeatFlux(timeSteps.size());
            Eigen::VectorXd probeCenter_gRec = probeCenter_trueHeatFlux;
            Eigen::VectorXd probeCenter_Ttrue = probeCenter_trueHeatFlux;
            Eigen::VectorXd probeCenter_Trec = probeCenter_trueHeatFlux;

            label cellprobeSide = mesh.findCell(probeSide);
            Eigen::VectorXd probeSide_trueHeatFlux(timeSteps.size());
            Eigen::VectorXd probeSide_gRec = probeSide_trueHeatFlux;

            label cellProbe2 = mesh.findCell(probe2);
            Eigen::VectorXd probe2_Ttrue(timeSteps.size());
            Eigen::VectorXd probe2_Trec = probe2_Ttrue;

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "Treconstructed",
                                      outputFolder);
            Info << "debug: Trec.size() = " << Trec.size() << endl;
            

            Info << "######### Postprocessing #########" <<
             endl;
            forAll(timeSteps, timeI)
            {
                volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                volScalarField gParametrizedField = list2Field(heatFlux[timeI]);
                ITHACAstream::exportSolution(trueHeatFluxField,
                                             std::to_string(timeSteps[timeI]), 
                                             outputFolder, "trueHeatFlux");
                ITHACAstream::exportSolution(Ttrue[timeI],
                        std::to_string(timeSteps[timeI]), outputFolder, "Ttrue");
                volScalarField relErrField = ITHACAutilities::relativeErrorField(
                        trueHeatFluxField, gParametrizedField);
                ITHACAstream::exportSolution(trueHeatFluxField,
                                             std::to_string(timeSteps[timeI]), 
                                             outputFolder, "trueHeatFlux");
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             outputFolder, "heatFluxRelErrField");
                relErrField = ITHACAutilities::relativeErrorField(Ttrue[timeI], 
                        Trec[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             outputFolder, "tempRelErrField");
                probeCenter_trueHeatFlux(timeI) = trueHeatFluxField[cellprobeCenter];
                probeCenter_gRec(timeI) = gParametrizedField[cellprobeCenter];
                probeCenter_Ttrue(timeI) = Ttrue[timeI][cellprobeCenter];
                probeCenter_Trec(timeI) = Trec[timeI][cellprobeCenter];
                probeSide_trueHeatFlux(timeI) = trueHeatFluxField[cellprobeSide];
                probeSide_gRec(timeI) = gParametrizedField[cellprobeSide];
                probe2_Ttrue(timeI) = Ttrue[timeI][cellProbe2];
                probe2_Trec(timeI) = Trec[timeI][cellProbe2];
                scalar EPS = 1e-6;
                heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, 
                        gParametrizedField, "hotSide");
                heatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                          gParametrizedField,
                                          "hotSide");
                trueHeatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, 
                        trueHeatFluxField, "hotSide");
                trueHeatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                          trueHeatFluxField,
                                          "hotSide");
                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            Info << "\nExporting" << endl;
            ITHACAstream::exportMatrix(heatFluxL2norm, "heatFluxL2norm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "heatFluxLinfNorm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(trueHeatFluxL2norm, "trueHeatFluxL2norm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(trueHeatFluxLinfNorm, "trueHeatFluxLinfNorm", 
                    "eigen", outputFolder);
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = 
                Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(
                    timeSteps_eigen, "timeSteps", "eigen", outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(probeCenter_trueHeatFlux, 
                    "probeCenter_trueHeatFlux", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_gRec, "probeCenter_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Ttrue, "probeCenter_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Trec, "probeCenter_Trec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_trueHeatFlux, "probeSide_trueHeatFlux", 
                    "eigen", outputFolder);
            ITHACAstream::exportMatrix(probeSide_gRec, "probeSide_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probe2_Ttrue, "probe2_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probe2_Trec, "probe2_Trec", "eigen", outputFolder);
            Info << "Ended Exporting" << endl;
        }

        void inverseProblemPostProcess(word outputFolder, word outputFolderFULL)
        {
            fvMesh& mesh(_mesh());
            Eigen::VectorXd heatFluxL2norm(timeSteps.size());
            Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());

            forAll(timeSteps, timeI)
            {
                volScalarField gParametrizedField = list2Field(heatFlux[timeI]);

                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(timeSteps[timeI]), 
                                             outputFolder, "gRecRED");
                heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, gParametrizedField,
                                        "hotSide");
                heatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                          gParametrizedField,
                                          "hotSide");
                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            ITHACAstream::exportMatrix(heatFluxL2norm, "heatFluxL2norm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "heatFluxLinfNorm", "eigen",
                                       outputFolder);
        }

        void noisePostProcess(word outputFolder, label testI)
        {
            word folder = "./ITHACAoutput/noiseTestVect/";
            fvMesh& mesh(_mesh());

            label cellprobeCenter = mesh.findCell(probeCenter);
            label cellprobeSide = mesh.findCell(probeSide);
            probeCenter_trueHeatFlux.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_gRec.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_Ttrue.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_Trec.conservativeResize(timeSteps.size(), testI + 1);

            probeSide_trueHeatFlux.conservativeResize(timeSteps.size(), testI + 1);
            probeSide_gRec.conservativeResize(timeSteps.size(), testI + 1);

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "Treconstructed",
                                      outputFolder);
            Info << "debug: Trec.size() = " << Trec.size() << endl;
            

            Info << "######### Postprocessing #########" <<
             endl;
            forAll(timeSteps, timeI)
            {
                volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                volScalarField gParametrizedField = list2Field(heatFlux[timeI]);
                volScalarField gDiffField = gParametrizedField - trueHeatFluxField;
                probeCenter_trueHeatFlux(timeI, testI) = trueHeatFluxField[cellprobeCenter];
                probeCenter_gRec(timeI, testI) = gParametrizedField[cellprobeCenter];
                probeCenter_Ttrue(timeI, testI) = Ttrue[timeI][cellprobeCenter];
                probeCenter_Trec(timeI, testI) = Trec[timeI][cellprobeCenter];
                probeSide_trueHeatFlux(timeI, testI) = trueHeatFluxField[cellprobeSide];
                probeSide_gRec(timeI, testI) = gParametrizedField[cellprobeSide];

                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = 
                Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(
                    timeSteps_eigen, "timeSteps", "eigen", outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(probeCenter_trueHeatFlux, 
                    "probeCenter_trueHeatFlux", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probeCenter_gRec, "probeCenter_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Ttrue, "probeCenter_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Trec, "probeCenter_Trec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_trueHeatFlux, "probeSide_trueHeatFlux", 
                    "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_gRec, "probeSide_gRec", "eigen", 
                    outputFolder);
        }

        //--------------------------------------------------------------------------
        /// Set parameter prior projecting initial gTrue on the parameterized space
        List<scalar> projectHeatFlux(List<scalar> _heatFlux)
        {
            Eigen::MatrixXd Temp(heatFluxSpaceBasis[0].size(),
                    heatFluxSpaceBasis.size());
            forAll(heatFluxSpaceBasis, baseI)
            {
                Temp.col(baseI) = Foam2Eigen::List2EigenMatrix(heatFluxSpaceBasis[baseI]);
            }

            Eigen::MatrixXd Btemp = Foam2Eigen::List2EigenMatrix(_heatFlux);
            Eigen::MatrixXd B = Temp.transpose() * Btemp;
            Temp = Temp.transpose() * Temp;
            Eigen::MatrixXd out = Temp.fullPivLu().solve(B);
            return Foam2Eigen::EigenMatrix2List(out);
        };

        void reconstructionTest(volScalarField _initialField, word _outputFolder)
        {
            Info << endl << endl << 
                "############ Testing T recostruction" << endl << endl;
            timeSampleI = 0;
            heatFluxWeightsOld = projectHeatFlux(trueHeatFlux[0]);
            recT.resize(0);
            recT.append(_initialField.clone());
            solT.resize(0);
            solT.append(_initialField.clone());


            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttime[NtimeStepsBetweenSamples -1]);
                     heatFluxWeightsOld = heatFluxWeights;
                }

                //Project true heat flux
                heatFluxWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);
                Info << "Weights = \n" << heatFluxWeights << endl;
                updateHeatFlux(heatFluxWeights);

                //Solve T0
                solveT0(_initialField);
                
                //ReconstructT
                reconstrucT(_outputFolder);
                forAll(Ttime, tI)
                {
                     recT.append(Ttime[tI].clone());
                }

                //Solve T
                solveT(_initialField, _outputFolder);
                forAll(Ttime, tI)
                {
                     solT.append(Ttime[tI].clone());
                }

                timeSampleI++;
            }

            //Compute relative error
            Info << "######### Postprocessing #########" <<
             endl;
            std::string message = "The recT.size() = " + std::to_string(recT.size())
                          + " while Ttrue.size is "
                          + std::to_string(Ttrue.size());

            M_Assert(recT.size() == Ttrue.size(), "Wrong size");
            M_Assert(solT.size() == Ttrue.size(), "Wrong size");
            forAll(timeSteps, timeI)
            {
                volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                volScalarField gParametrizedField = list2Field(heatFlux[timeI]);
                volScalarField relErrField = ITHACAutilities::relativeErrorField(
                        trueHeatFluxField, gParametrizedField);
                ITHACAstream::exportSolution(trueHeatFluxField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "trueHeatFlux");
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "heatFluxRelErrField");
                ITHACAstream::exportSolution(Ttrue[timeI],
                        std::to_string(timeSteps[timeI]), _outputFolder, "Ttrue");
                relErrField = ITHACAutilities::relativeErrorField(Ttrue[timeI], 
                        recT[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "tempRelErrField");
                relErrField = ITHACAutilities::relativeErrorField(Ttrue[timeI], 
                        solT[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "solTempRelErrField");
                relErrField = ITHACAutilities::relativeErrorField(solT[timeI], 
                        recT[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "solTrecTrelErrField");
                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            Info << "\nExporting" << endl;
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = 
                Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(
                    timeSteps_eigen, "timeSteps", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen",
                                       _outputFolder);
            Info << "Ended Exporting" << endl;
        }

        void pointsReconstructionTest(vector _point, volScalarField _initialField, 
                word _outputFolder)
        {
            Info << endl << endl << 
                "############ Testing poitwise T recostruction" << endl << endl;
            timeSampleI = 0;
            heatFluxWeightsOld = projectHeatFlux(trueHeatFlux[0]);
            recT.resize(0);
            recT.append(_initialField.clone());
            Eigen::VectorXd recT_vec(timeSamplesNum);
            Eigen::VectorXd toCheckT_vec(timeSamplesNum);
            Eigen::VectorXd Ttrue_vec(timeSamplesNum);
            Eigen::VectorXd recT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd toCheckT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd toCheckTrecT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd thetaT_vec(timeSamplesNum);

            Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
            Eigen::MatrixXd ThetaTilde = Theta + Temp;
            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttime[NtimeStepsBetweenSamples -1]);
                     heatFluxWeightsOld = heatFluxWeights;
                }

                //Project true heat flux
                heatFluxWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);
                Info << "Weights = \n" << heatFluxWeights << endl;
                updateHeatFlux(heatFluxWeights);

                //Solve T0
                solveT0(_initialField);
                
                //ReconstructT
                reconstrucT(_outputFolder);

                //for(int timeI = 0; timeI < NtimeStepsBetweenSamples; timeI++)
                //{
                //    label realTimeStep = timeI + 
                //        NtimeStepsBetweenSamples * timeSampleI;
                //    Info << "debug: realTimeStep = " << realTimeStep << endl;
                //    Ttrue_vec(realTimeStep) = fieldValueAtPoint(Ttrue[realTimeStep + 1], 
                //            _point);
                //    recT_vec(realTimeStep) = fieldValueAtPoint(Ttime[timeI], _point);
                //    toCheckT_vec(realTimeStep) = reconstrucT(_point, timeI);
                //}
                
                Ttrue_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttrue[samplingSteps[timeSampleI]], _point);
                recT_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttime[NtimeStepsBetweenSamples - 1], _point);
                toCheckT_vec(timeSampleI) = reconstrucT(_point);

                //Eigen::VectorXd w(heatFluxWeights.size());
                //forAll(heatFluxWeights, wI)
                //{
                //    w(wI) = heatFluxWeights[wI];
                //}
                //Eigen::VectorXd wOld(heatFluxWeights.size());
                //forAll(heatFluxWeights, wI)
                //{
                //    wOld(wI) = heatFluxWeightsOld[wI];
                //}
                //forAll(Ttime, tI)
                //{
                //     recT.append(Ttime[tI].clone());
                //     Eigen::VectorXd T0_vec = fieldValueAtThermocouples(T0_time[tI]);
                //     label realTimeI = samplingSteps[timeSampleI - 1] + tI;
                //     thetaT_vec.col(realTimeI) = ThetaTilde * w - addSol + T0_vec
                //         - Temp * wOld;
                //}

                timeSampleI++;
            }

            //Compute relative error
            Info << "######### Postprocessing #########" <<
             endl;
            std::string message = "The recT.size() = " + std::to_string(recT.size())
                          + " while Ttrue.size is "
                          + std::to_string(Ttrue.size());

            for(int timeI = 0; timeI < timeSamplesNum; timeI++)
            {
                recT_relErr_vec(timeI) = std::abs(Ttrue_vec(timeI) - recT_vec(timeI)) / 
                    Ttrue_vec(timeI);
                toCheckT_relErr_vec(timeI) = std::abs(Ttrue_vec(timeI) - 
                        toCheckT_vec(timeI)) / Ttrue_vec(timeI);
                toCheckTrecT_relErr_vec(timeI) = std::abs(recT_vec(timeI) - 
                        toCheckT_vec(timeI)) /  recT_vec(timeI);
                ITHACAstream::printProgress(double(timeI) / Ntimes);
            }
            Info << "\nExporting" << endl;
            ITHACAstream::exportMatrix(
                    recT_vec, "recT_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    Ttrue_vec, "Ttrue_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckT_vec, "toCheckT_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    recT_relErr_vec, "recT_relErr_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckT_relErr_vec, "toCheckT_relErr_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckTrecT_relErr_vec, "toCheckTrecT_relErr_vec", "eigen", 
                    _outputFolder);
            Info << "Ended Exporting" << endl;
        }

        void thermocoupleReconstructionTest(label _TCindex, volScalarField _initialField, 
                word _outputFolder)
        {
            Info << endl << endl << 
                "############ Testing poitwise T recostruction" << endl << endl;
            timeSampleI = 0;
            heatFluxWeightsOld = projectHeatFlux(trueHeatFlux[0]);
            recT.resize(0);
            recT.append(_initialField.clone());
            Eigen::VectorXd recT_vec(timeSamplesNum);
            Eigen::VectorXd toCheckT_vec(timeSamplesNum);
            Eigen::VectorXd Ttrue_vec(timeSamplesNum);
            Eigen::VectorXd recT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd toCheckT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd toCheckTrecT_relErr_vec(timeSamplesNum);
            Eigen::VectorXd thetaT_vec(timeSamplesNum);
            vector _point = thermocouplesPos[_TCindex];

            Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
            Eigen::MatrixXd ThetaTilde = Theta + Temp;
            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttime[NtimeStepsBetweenSamples -1]);
                     heatFluxWeightsOld = heatFluxWeights;
                }

                //Project true heat flux
                heatFluxWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);
                Info << "Weights = \n" << heatFluxWeights << endl;
                updateHeatFlux(heatFluxWeights);

                //Solve T0
                solveT0(_initialField);
                
                //ReconstructT
                reconstrucT(_outputFolder);
                
                Ttrue_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttrue[samplingSteps[timeSampleI]], _point);
                recT_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttime[NtimeStepsBetweenSamples - 1], _point);
                toCheckT_vec(timeSampleI) = reconstrucTatThermocouple(_TCindex);

                timeSampleI++;
            }

            //Compute relative error
            Info << "######### Postprocessing #########" <<
             endl;
            std::string message = "The recT.size() = " + std::to_string(recT.size())
                          + " while Ttrue.size is "
                          + std::to_string(Ttrue.size());

            for(int timeI = 0; timeI < timeSamplesNum; timeI++)
            {
                recT_relErr_vec(timeI) = std::abs(Ttrue_vec(timeI) - recT_vec(timeI)) / 
                    Ttrue_vec(timeI);
                toCheckT_relErr_vec(timeI) = std::abs(Ttrue_vec(timeI) - 
                        toCheckT_vec(timeI)) / Ttrue_vec(timeI);
                toCheckTrecT_relErr_vec(timeI) = std::abs(recT_vec(timeI) - 
                        toCheckT_vec(timeI)) /  recT_vec(timeI);
                ITHACAstream::printProgress(double(timeI) / Ntimes);
            }
            Info << "\nExporting" << endl;
            ITHACAstream::exportMatrix(
                    recT_vec, "recT_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    Ttrue_vec, "Ttrue_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckT_vec, "toCheckT_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    recT_relErr_vec, "recT_relErr_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckT_relErr_vec, "toCheckT_relErr_vec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    toCheckTrecT_relErr_vec, "toCheckTrecT_relErr_vec", "eigen", 
                    _outputFolder);
            Info << "Ended Exporting" << endl;
        }


        void linearSystemTest(volScalarField _initialField, label _Nimplicit, 
                word _outputFolder)
        {
            Info << endl << endl << 
                "############ Testing weights linear system" << endl << endl;
            timeSampleI = 0;

            Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
            Eigen::MatrixXd ThetaTilde = Theta + Temp;
            Eigen::MatrixXd residual(thermocouplesNum, timeSamplesNum);
            Eigen::VectorXd residualNorm(timeSamplesNum);
            Eigen::VectorXd weightsError(timeSamplesNum);
            List<scalar> trueWeights = 
                projectHeatFlux(trueHeatFlux[0]);
            heatFluxWeights = trueWeights;
            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                heatFluxWeightsOld = heatFluxWeights;
                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttrue[samplingSteps[timeSampleI - 1]]);
                     //ITHACAutilities::assignIF(_initialField, 
                     //        Ttime[NtimeStepsBetweenSamples -1]);
                }

                //Project true heat flux
                trueWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);
                Info << "True weights = \n" << trueWeights << endl; 
                Eigen::VectorXd trueWeights_Eig =
                    Foam2Eigen::List2EigenMatrix(trueWeights);

                //Solve T0
                solveT0(_initialField);

                Eigen::VectorXd heatFluxWeightsOld_Eig =
                    Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
                TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, 
                        thermocouplesNum *
                        NsamplesWindow);

                List<Eigen::MatrixXd> linSys;
                linSys.resize(2);
        
                Eigen::VectorXd weigths;
                if(linearBasis)
                {
                    linSys[0] = ThetaTilde.transpose() * ThetaTilde;
                    linSys[1] = ThetaTilde.transpose() * (TmeasShort + addSol +
                            Temp * heatFluxWeightsOld_Eig - T0_vector);
                    weigths = solveLinSys(linSys);
                }
                else
                {
                    linSys[0] = Theta.transpose() * Theta;
                    linSys[1] = Theta.transpose() * (TmeasShort + addSol - T0_vector);
                    weigths = solveLinSys(linSys);
                }

                Eigen::VectorXd tempError = (trueWeights_Eig - weigths).cwiseAbs();
                scalar EPS = 1e-36;
                for(int i = 0; i < tempError.size(); i++)
                {
                    scalar absTrue = std::abs(trueWeights_Eig(i));
                    if(absTrue < EPS)
                    {
                        tempError(i) =  tempError(i) / EPS;
                    }
                    else
                    {
                        tempError(i) =  tempError(i) / absTrue;
                    }
                }
                weightsError(timeSampleI) = tempError.norm();

                heatFluxWeights.resize(weigths.size());
                forAll(heatFluxWeights, weightI)
                {
                    heatFluxWeights[weightI] = weigths(weightI);
                }
                Info << "Computed weights = \n" << heatFluxWeights << endl; 
                Eigen::VectorXd R = ThetaTilde * weigths - addSol + T0_vector -
                    Temp * heatFluxWeightsOld_Eig - TmeasShort;
                residual.col(timeSampleI) = R;
                residualNorm(timeSampleI) = R.norm();

                updateHeatFlux(heatFluxWeights);

                //ReconstructT
                reconstrucT(_outputFolder);

                timeSampleI++;
            }

            //Compute relative error
            Info << "######### Postprocessing #########" << endl;

            PtrList<volScalarField> gRec;
            ITHACAstream::read_fields(gRec, "gReconstructed",
                                      _outputFolder);

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "Treconstructed",
                                      _outputFolder);

            Info << "\n*********************************************************" << endl;
            Info << "\nStart ONLY PROJECTION" << endl;
            timeSampleI = 0;

            trueWeights = 
                projectHeatFlux(trueHeatFlux[0]);
            ITHACAutilities::assignIF(_initialField, 
                    Ttrue[0]);
            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                heatFluxWeightsOld = trueWeights;
                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttrue[samplingSteps[timeSampleI - 1]]);
                }

                //Project true heat flux
                trueWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);

                //Solve T0
                solveT0(_initialField);
                heatFluxWeights = trueWeights;
                updateHeatFlux(heatFluxWeights);

                //ReconstructT
                reconstrucT(_outputFolder);

                timeSampleI++;
            }

            PtrList<volScalarField> gProj;
            ITHACAstream::read_fields(gProj, "gReconstructed",
                                      _outputFolder);

            PtrList<volScalarField> Tproj;
            ITHACAstream::read_fields(Tproj, "Treconstructed",
                                      _outputFolder);
            Info << "\nEnd ONLY PROJECTION" << endl;
            Info << "\n*********************************************************" << endl;


            Eigen::VectorXd heatFluxRelErr_L2(Ntimes);
            Eigen::VectorXd heatFluxRelErr_Linf(Ntimes);
            Eigen::VectorXd timeVec(Ntimes);
            Eigen::VectorXd probeHeat_rec(Ntimes);
            Eigen::VectorXd probeHeat_true(Ntimes);
            fvMesh& mesh = _mesh();
            label cellprobeCenter = mesh.findCell(probeCenter);
            forAll(timeSteps, timeI)
            {
                volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                ITHACAstream::exportSolution(trueHeatFluxField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "trueHeatFlux");
                volScalarField relErrField = ITHACAutilities::relativeErrorField(
                        trueHeatFluxField, gRec[timeI]);
                if(timeI > 0)
                {
                    probeHeat_true(timeI-1) = 
                        trueHeatFluxField[cellprobeCenter];
                    probeHeat_rec(timeI-1) = gRec[timeI][cellprobeCenter];
                    heatFluxRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(trueHeatFluxField, gRec[timeI]);
                    heatFluxRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(trueHeatFluxField, gRec[timeI]);
                    timeVec(timeI-1) = timeSteps[timeI];
                }
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "heatFluxRelErrField");
                ITHACAstream::exportSolution(gRec[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "computedHeatFlux");
                ITHACAstream::exportSolution(Trec[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "computedTemp");
                relErrField = ITHACAutilities::relativeErrorField(
                        Ttrue[timeI], Trec[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "tempRelErrField");
                relErrField = ITHACAutilities::relativeErrorField(
                        trueHeatFluxField, gProj[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "projHeatRelErr");
                relErrField = ITHACAutilities::relativeErrorField(
                        gProj[timeI], gRec[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "projRecHeatRelErr");
                relErrField = ITHACAutilities::relativeErrorField(
                        Tproj[timeI], Trec[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "projRecTempRelErr");
                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            Info << "\nEND" << endl;

            Info << "\nExporting" << endl;
            ITHACAstream::exportMatrix(
                    residual, "residual", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    residualNorm, "residualNorm", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    weightsError, "weightsError", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    heatFluxRelErr_L2, "heatFluxRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    heatFluxRelErr_Linf, "heatFluxRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    timeVec, "timeVec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    probeHeat_true, "probeHeat_true", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    probeHeat_rec, "probeHeat_rec", "eigen", _outputFolder);
            //ITHACAstream::exportMatrix(
            //        relErr, "relErr", "eigen", _outputFolder);
            //ITHACAstream::exportMatrix(
            //        compWeights, "compWeigths", "eigen", _outputFolder);
            //ITHACAstream::exportMatrix(
            //        projWeights, "projWeigths", "eigen", _outputFolder);
            Info << "Ended Exporting" << endl;
        }

        void updateHeatFlux2(Eigen::VectorXd _weights, word _outputFolder)
        {

            List<scalar> weights(_weights.size());
            forAll(weights, i)
            {
                weights[i] = _weights(i);
            }
            M_Assert(weights.size() == Nbasis,                                          
                     "weigths size different from basis functions size");                 
            volScalarField& T = _T();
            label firstTimeI = timeSampleI * NtimeStepsBetweenSamples;
            int lastTimeStep = firstTimeI + NtimeStepsBetweenSamples;                     
            updateTimeHeatFlux(weights);
            if(timeSampleI == 0)
            {
                volScalarField gParametrizedField = list2Field(heatFlux[0]);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(timeSteps[0]),
                                             _outputFolder,
                                             "gLinSys");
            }

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
                volScalarField gParametrizedField = list2Field(heatFlux[timeI+1]);
                ITHACAstream::exportSolution(gParametrizedField,
                                             std::to_string(timeSteps[timeI+1]),
                                             _outputFolder,
                                             "gLinSys");
            }   
        }

        void parameterizedBC(word outputFolder, volScalarField _initialField,
                List<scalar> _initialWeights)
        {
            Info << endl << "OVERWRITTEN Using quasilinearity of direct problem :" << endl;
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

        void computeLinearHeatWeights(word outputFolder,
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
                     ITHACAutilities::assignIF(_initialField, 
                             Ttime[NtimeStepsBetweenSamples -1]);
                     heatFluxWeightsOld = heatFluxWeights;
                     //heatFluxWeightsOld = 
                     //    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI - 1]]);
                }
                solveT0(_initialField);
        
                TmeasShort = Tmeas.segment(thermocouplesNum * timeSampleI, 
                        thermocouplesNum *
                        NsamplesWindow);
                Eigen::MatrixXd Temp = (1.0 / timeSamplesDeltaT) * Theta_tau;
                Eigen::MatrixXd M = Theta + Temp;
                Eigen::VectorXd heatFluxWeightsOld_Eig =
                    Foam2Eigen::List2EigenMatrix(heatFluxWeightsOld);
                linSys[0] = M.transpose() * M;
                linSys[1] = M.transpose() * (TmeasShort + addSol +
                        Temp * heatFluxWeightsOld_Eig - T0_vector);
        
                Eigen::VectorXd weigths = solveLinSys(linSys);
        
        
                heatFluxWeights.resize(weigths.size());
                forAll(heatFluxWeights, weightI)
                {
                    heatFluxWeights[weightI] = weigths(weightI);
                }
                Info << "Weights = \n" << heatFluxWeights << endl;
                updateHeatFlux(heatFluxWeights);
                label verbose = 0;
                parameterizedBC_postProcess(linSys, weigths, _initialField, outputFolder, 
                        verbose);
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

};

#endif

