class DanieliMold_unsteady: public inverseHeatTransferProblem
{
    public:
        explicit DanieliMold_unsteady(int argc, char* argv[])
            :
            inverseHeatTransferProblem(argc, argv)
        {
            fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
            timeSamplesNum = 10;
            timeSamplesT0 = 0.0;
            timeSamplesDeltaT = 1.0;
            samplingTime.resize(timeSamplesNum);
            forAll(samplingTime, timeI)
            {
                samplingTime[timeI] = timeSamplesT0 + timeI * timeSamplesDeltaT;
            }
            Info << "debug: samplingTimes = " << samplingTime << endl;
        }

        scalar timeSamplesDeltaT;
        scalar timeSamplesT0;

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

        volScalarField analyticalErrorField (volScalarField& field1,
                                             volScalarField& field2)
        {
            volScalarField fieldDiff(field1 - field2);
            return fieldDiff;
        }

        void readThermocouples()
        {
            if (!thermocouplesRead)
            {
                word fileName = "./thermocouplesCellsID";

                if (ITHACAutilities::check_file(fileName + "_mat.txt"))
                {
                    Info << "Reading thermocouples cells from file" << endl;
                    thermocouplesCellID = EigenMatrix2List(ITHACAstream::readMatrix(
                            fileName + "_mat.txt").cast <int> ());
                    thermocouplesNum = thermocouplesCellID.size();
                }
                else
                {
                    Info << "Defining positions of thermocouples" << endl;
                    fvMesh& mesh = _mesh();
                    volScalarField& T = _T();
                    Info << "Reading thermocouplesCoords file " << endl;
                    Eigen::MatrixXd thermocouplesCoords =
                        ITHACAstream::readMatrix("./thermocouplesCoords.txt");
                    thermocouplesCellID.resize(thermocouplesCoords.rows());
                    thermocouplesNum = thermocouplesCellID.size();

                    for (label tcI = 0; tcI < thermocouplesCoords.rows(); tcI++)
                    {
                        point tcCoord(thermocouplesCoords(tcI, 0), thermocouplesCoords(tcI, 1),
                                      thermocouplesCoords(tcI, 2));
                        thermocouplesCellID[tcI] = mesh.findCell(tcCoord);
                    }

                    volScalarField thermocouplesField(T);
                    ITHACAutilities::assignIF(thermocouplesField, homogeneousBC);
                    forAll(thermocouplesCellID, tcI)
                    {
                        thermocouplesField.ref()[thermocouplesCellID[tcI]] = 1;
                    }
                    ITHACAstream::exportSolution(thermocouplesField, "1", "./ITHACAoutput/debug/",
                                                 "thermocouplesField,");
                    Eigen::MatrixXi thermocouplesCellID_eigen = List2EigenMatrix(
                                thermocouplesCellID);
                    ITHACAstream::exportMatrix(thermocouplesCellID_eigen, fileName,
                                               "eigen", "./");
                    ////Define thermocouples plane
                    //defineThermocouplesPlane();
                    //Info << "Exporting debug field" << endl;
                    //ITHACAutilities::assignIF(T, homogeneousBC);
                    //forAll(thermocouplesPlaneCellID, cellI)
                    //{
                    //    T.ref()[thermocouplesPlaneCellID[cellI]] = 1000;
                    //}
                    //ITHACAstream::exportSolution(T, "1", "./ITHACAoutput/debug/",
                    //                             "thermocouplesSurface");
                }

                thermocouplesRead = 1;
            }
            else
            {
                WarningInFunction << "readThermocouples function called twice." << endl;
                WarningInFunction << "I am not doing the second reading." << endl;
            }
        }

        void assignBC()
        {
            //  I must rewrite this function considering that gTrue is now a Lis of List
            //    fvMesh& mesh = _mesh();
            //    volScalarField& T(_T());
            //    // Setting BC at the cold side
            //    coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            //    label coldSideSize = T.boundaryField()[coldSide_ind].size();
            //    Tf.resize(coldSideSize);
            //    set_valueFraction();
            //    forAll(Tf, faceI)
            //    {
            //        scalar faceZ =
            //            mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
            //        Tf[faceI] = 320 - faceZ * 9;
            //    }
            //    // Setting BC at hotSide
            //    hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            //    label hotSideSize = T.boundaryField()[hotSide_ind].size();
            //    gTrue.resize(hotSideSize);
            //    forAll(gTrue, faceI)
            //    {
            //        scalar faceX =
            //            mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
            //        scalar faceZ =
            //            mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
            //        scalar faceArea =
            //      mesh.magSf().boundaryField()[hotSide_ind][faceI];
            //        gTrue[faceI] = - 1e6;
            //        //gTrue[faceI] = 0;
            //    }
            //    set_valueFraction();
            //    List<scalar> RobinBC = Tf;
            //    forAll(mesh.boundaryMesh(), patchI)
            //    {
            //        if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
            //        {
            //            ITHACAutilities::assignMixedBC(T, patchI, RobinBC, refGrad, valueFraction);
            //        }
            //        else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
            //        {
            //            ITHACAutilities::assignBC(T, patchI, - gTrue / k);
            //        }
            //        else
            //        {
            //            ITHACAutilities::assignBC(T, patchI, homogeneousBC);
            //        }
            //    }
            //
        }
        void assignT0_IF(volScalarField& T0_field)
        {
        }

        void updateBC()
        {
            //  I must rewrite this function considering that gTrue is now a Lis of List
            // fvMesh& mesh = _mesh();
            // volScalarField& T(_T());
            // // Setting BC at hotSide
            // hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            // label hotSideSize = T.boundaryField()[hotSide_ind].size();
            // gTrue.resize(hotSideSize);
            // forAll(gTrue, faceI)
            // {
            //     scalar faceX =
            //         mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
            //     scalar faceZ =
            //         mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
            //     gTrue[faceI] *= 0.1;
            // }
            // set_valueFraction();
            // List<scalar> RobinBC = Tf;
            // forAll(mesh.boundaryMesh(), patchI)
            // {
            //     if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
            //     {
            //         ITHACAutilities::assignBC(T, patchI, - gTrue / k);
            //     }
            // }
        }

        //void solveSteady()
        //{
        //    //Reinitializing
        //    restart();

        //    assignBC();
        //    //updateBC();
        //    volScalarField& T = _T();
        //    simpleControl& simple = _simple();

        //    while (simple.loop())
        //    {
        //        while (simple.correctNonOrthogonal())
        //        {
        //            fvScalarMatrix TEqn
        //            (
        //                fvm::laplacian(DT, T)
        //            );
        //            TEqn.solve();
        //        }
        //    }
        //}

        //void steadyTest( )
        //{
        //    word folderOffline = "./ITHACAoutput/test/";
        //    label savingTime = 1;
        //    solveSteady();
        //    volScalarField& T(_T());
        //    ITHACAstream::exportSolution(T, std::to_string(savingTime),
        //                     folderOffline,
        //                     "Tsteady");
        //    //timeStep = 0;
        //    //List<scalar> temperatureOverTime;
        //    //for(label i = 0; i < heatFlux.rows(); i++)
        //    //{
        //    //    Foam::Time& runTime = _runTime();
        //    //    solveSteady();
        //    //    volScalarField& T(_T());
        //    //    ITHACAstream::exportSolution(T, std::to_string(savingTime),
        //    //                     folderOffline,
        //    //                     "Tsteady");
        //    //    temperatureOverTime.append(T.internalField()[thermocouplesCellID[0]]);
        //    //    savingTime++;
        //    //    timeStep++;
        //    //}
        //    //Eigen::MatrixXd temperatureOverTime_eigen = List2EigenMatrix(temperatureOverTime);
        //    //ITHACAstream::exportMatrix(temperatureOverTime_eigen, "temperatureOverTime_steady", "eigen", folderOffline);
        //    Info<< "End steady\n" << endl;
        //}

        void solveUnsteady()
        {
            restart();
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());
            Foam::Time& runTime = _runTime();
            const char* folderOffline = "./ITHACAoutput/test/";
            Info << "Thermal diffusivity = " << diffusivity << " m2/s" << endl;
            simpleControl& simple = _simple();
            fv::options& fvOptions(_fvOptions());
            assignBC();
            bool updatedBC = 0;
            solve("direct", folderOffline);
            Info << "End\n" << endl;
        }
};
