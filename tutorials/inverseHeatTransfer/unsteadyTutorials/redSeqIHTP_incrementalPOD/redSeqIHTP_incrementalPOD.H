#ifndef redSeqIHTP_incrementalPOD_H
#define redSeqIHTP_incrementalPOD_H

#include "ReducedSequentialIHTP_incrementalPOD.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"


class redSeqIHTP_incrementalPOD: public reducedSequentialIHTP_incrementalPOD
{
    public:
        explicit redSeqIHTP_incrementalPOD(int argc, char* argv[])
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
        volScalarField initialField;

        Foam::vector probeCenter = Foam::vector(1, 0.0, 0.6);
        Foam::vector probeSide = Foam::vector(0.8, 0.0, 1.1);
        Foam::vector probe2 = Foam::vector(1.0, 0.02, 0.6);

        Eigen::MatrixXd probeCenter_gTrue;
        Eigen::MatrixXd probeCenter_gRec;
        Eigen::MatrixXd probeCenter_Ttrue;
        Eigen::MatrixXd probeCenter_Trec;

        Eigen::MatrixXd probeSide_gTrue;
        Eigen::MatrixXd probeSide_gRec;

        List<scalar> heatFlux_steady;
        scalar maxFrequency; // Hz

        List<int> EigenMatrix2List ( Eigen::MatrixXi matrix )
        {
            if (matrix.cols() == 1)
            {
                List<int> list(matrix.rows());

                for (label i = 0; i < matrix.rows(); i++)
                {
                    list[i] = matrix(i, 0);
                }

                return list;
            }
            else
            {
                Info << "I still must code this function, exiting" << endl;
                exit(11);
            }
        }

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


        void assignT0_IF(volScalarField& T0_field)
        {
            assignTrueIF();
            volScalarField& T(_T());
            ITHACAutilities::assignIF(T0_field, T);
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
                heatFlux_steady[faceI] = - k * (a * faceX * faceX + b * faceZ + c);
            }
        }

        void set_gTrue()
        {
            fvMesh& mesh(_mesh());
            volScalarField& T(_T());
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
            gTrue.resize(timeSteps.size());
            set_steadyHeatFlux();
            forAll(timeSteps, timeI)
            {
                // Setting gTrue
                gTrue[timeI].resize(hotSideSize);
                forAll(gTrue[timeI], faceI)
                {
                    scalar faceX =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                    scalar faceZ =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();

                    scalar g1 = - k * (b * faceZ * faceZ  + c);
                    scalar g2 = - k * c / (1 + (faceX - 1) * (faceX - 1) + faceZ * faceZ);
                    gTrue[timeI][faceI] = g1 + 
                        g1 / 2 * (std::sin ( 2 * M_PI * maxFrequency * 
                                    timeSteps[timeI] * timeSteps[timeI] / endTime)) +
                        g2 * std::exp(- 0.1 * timeSteps[timeI]);


                    //gTrue[timeI][faceI] = heatFlux_steady[faceI] + 
                    //    0.1 * heatFlux_steady[faceI] * (
                    //            std::sin ( 2 * M_PI * maxFrequency * timeSteps[timeI]));

                    //scalar g1 = - k * (a * faceX  + b * faceZ * faceZ  + c);
                    //gTrue[timeI][faceI] += g1 * (1 + 0.5 * (std::sin (
                    //                          2 * M_PI * (maxFrequency / endTime * 
                    //                              timeSteps[timeI]) * timeSteps[timeI])));

                    //if(timeSteps[timeI] < 10.0)
                    //{
                    //    scalar g1 = - k * (a * faceX  + b * faceZ * faceZ  + c);
                    //    gTrue[timeI][faceI] = g1 * (1 + 0.5 * (std::sin (
                    //                              2 * M_PI * (maxFrequency / endTime * timeSteps[timeI]) * timeSteps[timeI])));
                    //}
                    
                    //scalar g1 = - k * (a * faceX * faceX * std::sin (2 * M_PI * timeSteps[timeI]) + b * faceZ * std::sin (M_PI * timeSteps[timeI]) + c);
                    //gTrue[timeI][faceI] = g1 * (1 + 0.5 * (std::sin (
                    //                          2 * M_PI * (maxFrequency / endTime * timeSteps[timeI]) * timeSteps[timeI])));
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
                    ITHACAutilities::assignBC(T, patchI, - g[timeI] / k);
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
            g = gTrue;
            assignDirectBC(timeI);
            ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI]),
                                         directOutputFolder,
                                         "Ttrue");
            volScalarField gTrueField = list2Field(gTrue[timeI]);
            ITHACAstream::exportSolution(gTrueField,
                                         std::to_string(timeSteps[timeI]), directOutputFolder,
                                         "gTrue");
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
                gTrueField = list2Field(gTrue[timeI]);
                ITHACAstream::exportSolution(gTrueField,
                                             std::to_string(timeSteps[timeI]), 
                                             directOutputFolder,
                                             "gTrue");
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
            Info << "debug : g.size() = " << g.size() << endl;
            Info << "debug : gTrue.size() = " << gTrue.size() << endl;
            Info << "debug : Ttrue.size() = " << Ttrue.size() << endl;
            Info << "debug : timeSteps.size() = " << timeSteps.size() << endl;

            label cellprobeCenter = mesh.findCell(probeCenter);
            Eigen::VectorXd probeCenter_gTrue(timeSteps.size());
            Eigen::VectorXd probeCenter_gRec = probeCenter_gTrue;
            Eigen::VectorXd probeCenter_Ttrue = probeCenter_gTrue;
            Eigen::VectorXd probeCenter_Trec = probeCenter_gTrue;

            label cellprobeSide = mesh.findCell(probeSide);
            Eigen::VectorXd probeSide_gTrue(timeSteps.size());
            Eigen::VectorXd probeSide_gRec = probeSide_gTrue;

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
                volScalarField gTrueField = list2Field(gTrue[timeI]);
                volScalarField gParametrizedField = list2Field(g[timeI]);
                volScalarField gDiffField = gParametrizedField - gTrueField;
                ITHACAstream::exportSolution(gTrueField,
                                             std::to_string(timeSteps[timeI]), outputFolder,
                                             "gTrue");
                ITHACAstream::exportSolution(Ttrue[timeI],
                                             std::to_string(timeSteps[timeI]), outputFolder,
                                             "Ttrue");
                ITHACAstream::exportSolution(gDiffField,
                                             std::to_string(timeSteps[timeI]), outputFolder,
                                             "gDiffField");
                probeCenter_gTrue(timeI) = gTrueField[cellprobeCenter];
                probeCenter_gRec(timeI) = gParametrizedField[cellprobeCenter];
                probeCenter_Ttrue(timeI) = Ttrue[timeI][cellprobeCenter];
                probeCenter_Trec(timeI) = Trec[timeI][cellprobeCenter];
                probeSide_gTrue(timeI) = gTrueField[cellprobeSide];
                probeSide_gRec(timeI) = gParametrizedField[cellprobeSide];
                probe2_Ttrue(timeI) = Ttrue[timeI][cellProbe2];
                probe2_Trec(timeI) = Trec[timeI][cellProbe2];
                scalar EPS = 1e-6;
                //volScalarField relativeErrorField(gDiffField);

                //for (label i = 0; i < relativeErrorField.internalField().size(); i++)
                //{
                //    if (std::abs(gTrueField.ref()[i]) < EPS)
                //    {
                //        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
                //    }
                //    else
                //    {
                //        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) /
                //                                      std::abs(gTrueField.ref()[i]);
                //    }
                //}

                //ITHACAstream::exportSolution(relativeErrorField,
                //                             std::to_string(timeSteps[timeI]), outputFolder,
                //                             "relativeErrorField");
                heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, gParametrizedField,
                                        "hotSide");
                heatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                          gParametrizedField,
                                          "hotSide");
                trueHeatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, gTrueField,
                                       "hotSide");
                trueHeatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                          gTrueField,
                                          "hotSide");
                ITHACAstream::printProgress(double(timeI + 1) / timeSteps.size());
            }
            //for(int sI = 0; sI < timeSamplesNum; sI++)
            //{
            //    double avgL2 = 0;
            //    double avgLinf = 0;
            //    for(int i = 0; i < NtimeStepsBetweenSamples; i++)
            //    {
            //        int timeI = sI * NtimeStepsBetweenSamples + i;
            //        avgL2 += heatFluxL2norm(timeI);
            //        avgLinf += heatFluxLinfNorm(timeI);
            //    }
            //    avgL2 = avgL2 / NtimeStepsBetweenSamples;
            //    avgLinf = avgLinf / NtimeStepsBetweenSamples;
            //    for(int i = 0; i < NtimeStepsBetweenSamples; i++)
            //    {
            //        int timeI = sI * NtimeStepsBetweenSamples + i;
            //        avgRelErrL2norm(timeI) = avgL2;
            //        avgRelErrLinfNorm(timeI) = avgLinf;
            //    }
            //    ITHACAstream::printProgress(double(sI + 1) / timeSamplesNum);
            //}
            Info << "\nExporting" << endl;
            ITHACAstream::exportMatrix(heatFluxL2norm, "heatFluxL2norm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "heatFluxLinfNorm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(trueHeatFluxL2norm, "trueHeatFluxL2norm", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(trueHeatFluxLinfNorm, "trueHeatFluxLinfNorm", "eigen",
                                       outputFolder);
            //ITHACAstream::exportMatrix(heatFluxL2norm, "relErrL2norm", "eigen",
            //                           outputFolder);
            //ITHACAstream::exportMatrix(heatFluxLinfNorm, "relErrLinfNorm", "eigen",
            //                           outputFolder);
            //ITHACAstream::exportMatrix(avgRelErrL2norm, "avgRelErrL2norm", "eigen",
            //                           outputFolder);
            //ITHACAstream::exportMatrix(avgRelErrLinfNorm, "avgRelErrLinfNorm", "eigen",
            //                           outputFolder);
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = 
                Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(
                    timeSteps_eigen, "timeSteps", "eigen", outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen",
                                       outputFolder);
            ITHACAstream::exportMatrix(probeCenter_gTrue, "probeCenter_gTrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_gRec, "probeCenter_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Ttrue, "probeCenter_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Trec, "probeCenter_Trec", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probeSide_gTrue, "probeSide_gTrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_gRec, "probeSide_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probe2_Ttrue, "probe2_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probe2_Trec, "probe2_Trec", "eigen", outputFolder);
            Info << "\n Ended Exporting" << endl;
        }

        void inverseProblemPostProcess(word outputFolder, word outputFolderFULL)
        {
            fvMesh& mesh(_mesh());
            Eigen::VectorXd heatFluxL2norm(timeSteps.size());
            Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());

            forAll(timeSteps, timeI)
            {
                volScalarField gParametrizedField = list2Field(g[timeI]);

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
            probeCenter_gTrue.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_gRec.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_Ttrue.conservativeResize(timeSteps.size(), testI + 1);
            probeCenter_Trec.conservativeResize(timeSteps.size(), testI + 1);

            probeSide_gTrue.conservativeResize(timeSteps.size(), testI + 1);
            probeSide_gRec.conservativeResize(timeSteps.size(), testI + 1);

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "Treconstructed",
                                      outputFolder);
            Info << "debug: Trec.size() = " << Trec.size() << endl;
            

            Info << "######### Postprocessing #########" <<
             endl;
            forAll(timeSteps, timeI)
            {
                volScalarField gTrueField = list2Field(gTrue[timeI]);
                volScalarField gParametrizedField = list2Field(g[timeI]);
                volScalarField gDiffField = gParametrizedField - gTrueField;
                probeCenter_gTrue(timeI, testI) = gTrueField[cellprobeCenter];
                probeCenter_gRec(timeI, testI) = gParametrizedField[cellprobeCenter];
                probeCenter_Ttrue(timeI, testI) = Ttrue[timeI][cellprobeCenter];
                probeCenter_Trec(timeI, testI) = Trec[timeI][cellprobeCenter];
                probeSide_gTrue(timeI, testI) = gTrueField[cellprobeSide];
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
            ITHACAstream::exportMatrix(probeCenter_gTrue, "probeCenter_gTrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_gRec, "probeCenter_gRec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Ttrue, "probeCenter_Ttrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeCenter_Trec, "probeCenter_Trec", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_gTrue, "probeSide_gTrue", "eigen", 
                    outputFolder);
            ITHACAstream::exportMatrix(probeSide_gRec, "probeSide_gRec", "eigen", 
                    outputFolder);
        }
};

#endif

