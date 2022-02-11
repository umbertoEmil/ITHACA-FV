#ifndef constantBasisTest_ROM_H
#define constantBasisTest_ROM_H

#include "ReducedSequentialIHTP_incrementalPOD_constant.H"
#include "MUQ/Modeling/Distributions/Gaussian.h"


class constantBasisTest_ROM: public reducedSequentialIHTP_incrementalPOD_constant
{
    public:
        explicit constantBasisTest_ROM(int argc, char* argv[])
            :
            reducedSequentialIHTP_incrementalPOD_constant(argc, argv),
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
        bool linearTrueHeatFlux = 0;

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


                    if(linearTrueHeatFlux)
                    {
                        scalar g1 = - thermalCond * (b * faceZ * faceZ  + c);
                        trueHeatFlux[timeI][faceI] = timeGrad * timeSteps[timeI] * g1 + g1;
                    }
                    else
                    {
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
        /// 
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
                heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, 
                        gParametrizedField, "hotSide");
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
};

#endif

