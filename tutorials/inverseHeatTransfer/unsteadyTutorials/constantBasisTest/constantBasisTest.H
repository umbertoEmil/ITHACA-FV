#ifndef constantBasisTest_H
#define constantBasisTest_H

//#include "ReducedSequentialIHTP_incrementalPOD.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"

class constantBasisTest: public sequentialIHTP_constant
{
    public:
        explicit constantBasisTest(int argc, char* argv[])
            :
            sequentialIHTP_constant(argc, argv),
            initialField(_T())
        {
            fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
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

        //--------------------------------------------------------------------------
        /// 
        Eigen::MatrixXi List2EigenMatrix ( List<int> list )
        {
            Eigen::MatrixXi matrix(list.size(), 1);

            for (label i = 0; i < matrix.rows(); i++)
            {
                matrix(i, 0) = list[i];
            }

            return matrix;
        }

        //--------------------------------------------------------------------------
        /// 
        Eigen::MatrixXd List2EigenMatrix ( List<scalar> list )
        {
            Eigen::MatrixXd matrix(list.size(), 1);

            for (label i = 0; i < matrix.rows(); i++)
            {
                matrix(i, 0) = list[i];
            }

            return matrix;
        }

        //--------------------------------------------------------------------------
        /// 
        void assignTrueIF()
        {
            volScalarField& T(_T());

            for (label i = 0; i < T.internalField().size(); i++)
            {
                T.ref()[i] =  initialField.internalField()[i];
            }
        }

        //--------------------------------------------------------------------------
        /// 
        /// Set steady component of heat flux, g
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

        //--------------------------------------------------------------------------
        /// 
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

                    scalar g1 = - thermalCond * (b * faceZ * faceZ  + c);
                    trueHeatFlux[timeI][faceI] = timeGrad * timeSteps[timeI] * g1 + g1; 

                    //scalar g1 = - thermalCond * (b * faceZ * faceZ  + c);
                    //scalar g2 = - thermalCond * c / (1 + (faceX - 1) * (faceX - 1) +
                    //        faceZ * faceZ);
                    //trueHeatFlux[timeI][faceI] = g1 +
                    //    g1 / 2 * (std::sin ( 2 * M_PI * maxFrequency *
                    //                timeSteps[timeI] * timeSteps[timeI] / endTime)) +
                    //    g2 * std::exp(- 0.1 * timeSteps[timeI]);

                }
            }
        }

        //--------------------------------------------------------------------------
        /// 
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

        //--------------------------------------------------------------------------
        /// 
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

        //--------------------------------------------------------------------------
        /// 
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

        //--------------------------------------------------------------------------
        /// 
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

        //--------------------------------------------------------------------------
        /// 
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
            ITHACAstream::read_fields(Trec, reconstructedTemperatureName,
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

        //--------------------------------------------------------------------------
        /// 
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

                if(timeSampleI == 0)
                {
                    heatFluxWeightsGrad_old = heatFluxWeights * 0.0;
                }
                else
                {
                    heatFluxWeightsGrad_old = heatFluxWeightsGrad;
                }
                heatFluxWeightsGrad = (heatFluxWeights - heatFluxWeightsOld) / 
                    timeSamplesDeltaT;

                updateHeatFlux(heatFluxWeights);

                solveT_ic(_initialField);
                
                //ReconstructT
                reconstructT(_outputFolder);

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
            
            Eigen::VectorXd heatFluxRelErr_L2(Ntimes);
            Eigen::VectorXd heatFluxRelErr_Linf(Ntimes);
            Eigen::VectorXd TsolRelErr_L2(Ntimes);
            Eigen::VectorXd TsolRelErr_Linf(Ntimes);
            Eigen::VectorXd TrecRelErr_L2(Ntimes);
            Eigen::VectorXd TrecRelErr_Linf(Ntimes);
            Eigen::VectorXd TrecTsolRelErr_L2(Ntimes);
            Eigen::VectorXd TrecTsolRelErr_Linf(Ntimes);
            Eigen::VectorXd timeVec(Ntimes);
            forAll(timeSteps, timeI)
            {
                volScalarField trueHeatFluxField = list2Field(trueHeatFlux[timeI]);
                volScalarField gParametrizedField = list2Field(heatFlux[timeI]);
                volScalarField relErrField = ITHACAutilities::relativeErrorField(
                        trueHeatFluxField, gParametrizedField);
                if(timeI > 0)
                {
                    heatFluxRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(trueHeatFluxField, gParametrizedField);
                    heatFluxRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(trueHeatFluxField, 
                                gParametrizedField);
                    TsolRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(Ttrue[timeI], solT[timeI]);
                    TsolRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(Ttrue[timeI], solT[timeI]);  
                    TrecRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(Ttrue[timeI], recT[timeI]);
                    TrecRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(Ttrue[timeI], recT[timeI]);  
                    TrecTsolRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(solT[timeI], recT[timeI]);
                    TrecTsolRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(solT[timeI], recT[timeI]);  
                    timeVec(timeI-1) = timeSteps[timeI];
                }
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
            ITHACAstream::exportMatrix(
                    heatFluxRelErr_L2, "heatFluxRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    heatFluxRelErr_Linf, "heatFluxRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TsolRelErr_L2, "TsolRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TsolRelErr_Linf, "TsolRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TrecRelErr_L2, "TrecRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TrecRelErr_Linf, "TrecRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TrecTsolRelErr_L2, "TrecTsolRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    TrecTsolRelErr_Linf, "TrecTsolRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    timeVec, "timeVec", "eigen", _outputFolder);
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = 
                Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(
                    timeSteps_eigen, "timeSteps", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen",
                                       _outputFolder);
            Info << "Ended Exporting" << endl;
        }

        //--------------------------------------------------------------------------
        /// 
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

                if(timeSampleI == 0)
                {
                    heatFluxWeightsGrad_old = heatFluxWeights * 0.0;
                }
                else
                {
                    heatFluxWeightsGrad_old = heatFluxWeightsGrad;
                }
                heatFluxWeightsGrad = (heatFluxWeights - heatFluxWeightsOld) / 
                    timeSamplesDeltaT;

                updateHeatFlux(heatFluxWeights);

                solveT_ic(_initialField);
                
                //ReconstructT
                reconstructT(_outputFolder);

                T_ic_ready = 1;

                toCheckT_vec(timeSampleI) = reconstructTatThermocouple(_TCindex);
                
                Ttrue_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttrue[samplingSteps[timeSampleI]], _point);
                recT_vec(timeSampleI) = 
                    fieldValueAtPoint(Ttime[NtimeStepsBetweenSamples - 1], _point);

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

        //--------------------------------------------------------------------------
        /// 
        void inverseProblemTest(volScalarField _initialField, List<scalar> _initialWeights,
                word _outputFolder)
        {
            Info << endl << endl << 
                "############ Inverse problem test ############" << endl << endl;

            computeHeatFluxWeights(_initialField, _initialWeights, _outputFolder);

            //Compute relative error
            Info << "######### Postprocessing #########" << endl;

            PtrList<volScalarField> gRec;
            ITHACAstream::read_fields(gRec, reconstructedHeatFluxName,
                                      _outputFolder);

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, reconstructedTemperatureName,
                                      _outputFolder);

            Info << "\n*********************************************************" << endl;
            Info << "\nStart ONLY PROJECTION" << endl;
            timeSampleI = 0;

            ITHACAutilities::assignIF(_initialField, 
                    Ttrue[0]);
            heatFluxWeights = _initialWeights;
            List<scalar> trueWeights = _initialWeights;
            while(timeSampleI < timeSamplesNum)
            {
                Info << "\nTime sample " << timeSampleI + 1 << endl;

                if(timeSampleI > 0)
                {
                     ITHACAutilities::assignIF(_initialField, 
                             Ttrue[samplingSteps[timeSampleI - 1]]);
                }

                heatFluxWeightsOld = trueWeights;

                //Project true heat flux
                trueWeights = 
                    projectHeatFlux(trueHeatFlux[samplingSteps[timeSampleI]]);
                Eigen::VectorXd trueWeights_Eig =
                    Foam2Eigen::List2EigenMatrix(trueWeights);

                solveT_ic(_initialField);

                heatFluxWeights = trueWeights;
                updateHeatFlux(heatFluxWeights);

                //ReconstructT
                reconstructT(_outputFolder);

                timeSampleI++;
            }

            PtrList<volScalarField> gProj;
            ITHACAstream::read_fields(gProj, reconstructedHeatFluxName,
                                      _outputFolder);

            PtrList<volScalarField> Tproj;
            ITHACAstream::read_fields(Tproj, reconstructedTemperatureName,
                                      _outputFolder);
            Info << "\nEnd ONLY PROJECTION" << endl;
            Info << "\n*********************************************************" << endl;


            Eigen::VectorXd heatFluxRelErr_L2(Ntimes);
            Eigen::VectorXd heatFluxRelErr_Linf(Ntimes);
            Eigen::VectorXd projHeatFluxRelErr_L2(Ntimes);
            Eigen::VectorXd projHeatFluxRelErr_Linf(Ntimes);
            Eigen::VectorXd projRecHeatFluxRelErr_L2(Ntimes);
            Eigen::VectorXd projRecHeatFluxRelErr_Linf(Ntimes);
            Eigen::VectorXd timeVec(Ntimes);
            Eigen::VectorXd probeHeat_rec(Ntimes);
            Eigen::VectorXd probeHeat_proj(Ntimes);
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
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "heatFluxRelErrField");
                if(timeI > 0)
                {
                    probeHeat_true(timeI-1) = 
                        trueHeatFluxField[cellprobeCenter];
                    probeHeat_rec(timeI-1) = gRec[timeI][cellprobeCenter];
                    probeHeat_proj(timeI-1) = gProj[timeI][cellprobeCenter];
                    heatFluxRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(trueHeatFluxField, gRec[timeI]);
                    heatFluxRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(trueHeatFluxField, gRec[timeI]);
                    projHeatFluxRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(trueHeatFluxField, gProj[timeI]);
                    projHeatFluxRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(trueHeatFluxField, gProj[timeI]);
                    projRecHeatFluxRelErr_L2(timeI-1) = 
                        ITHACAutilities::errorL2Rel(gProj[timeI], gRec[timeI]);
                    projRecHeatFluxRelErr_Linf(timeI-1) = 
                        ITHACAutilities::errorLinfRel(gProj[timeI], gRec[timeI]);
                    timeVec(timeI-1) = timeSteps[timeI];
                }
                ITHACAstream::exportSolution(gRec[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, reconstructedHeatFluxName);
                ITHACAstream::exportSolution(Trec[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, reconstructedTemperatureName);
                ITHACAstream::exportSolution(gProj[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "projHeatFlux");
                ITHACAstream::exportSolution(Tproj[timeI],
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "projTemp");
                relErrField = ITHACAutilities::relativeErrorField(
                        Ttrue[timeI], Trec[timeI]);
                ITHACAstream::exportSolution(relErrField,
                                             std::to_string(timeSteps[timeI]), 
                                             _outputFolder, "recTempRelErrField");
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
                    heatFluxRelErr_L2, "heatFluxRelErr_L2", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    heatFluxRelErr_Linf, "heatFluxRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(projHeatFluxRelErr_L2, "projHeatFluxRelErr_L2", 
                    "eigen", _outputFolder);
            ITHACAstream::exportMatrix( projHeatFluxRelErr_Linf, 
                    "projHeatFluxRelErr_Linf", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    timeVec, "timeVec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    probeHeat_true, "probeHeat_true", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    probeHeat_rec, "probeHeat_rec", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    probeHeat_proj, "probeHeat_proj", "eigen", _outputFolder);
            ITHACAstream::exportMatrix(
                    projRecHeatFluxRelErr_L2, "projRecHeatFluxRelErr_L2", "eigen", 
                    _outputFolder);
            ITHACAstream::exportMatrix(
                    projRecHeatFluxRelErr_Linf, "projRecHeatFluxRelErr_Linf", "eigen", 
                    _outputFolder);
            Info << "Ended Exporting" << endl;
        }

};

#endif

