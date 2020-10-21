#ifndef numericalChirpTest_unsteady_H
#define numericalChirpTest_unsteady_H

class numericalChirpTest_unsteady: public inverseHeatTransferProblem
{
    public:
        explicit numericalChirpTest_unsteady(int argc, char* argv[])
            :
            inverseHeatTransferProblem(argc, argv),
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

	List<scalar> heatFlux_steady;
	scalar maxFrequency; // Hz

        List<int> EigenMatrix2List ( Eigen::MatrixXi matrix )
	{
	    if(matrix.cols() == 1)
	    {
	        List<int> list(matrix.rows());
	        for(label i = 0; i < matrix.rows(); i++)
	        {
	            list[i] = matrix(i,0);
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
	    Eigen::MatrixXi matrix(list.size(),1);
	    for(label i = 0; i < matrix.rows(); i++)
	    {
	        matrix(i,0) = list[i];
	    }
	    return matrix;
	}

        Eigen::MatrixXd List2EigenMatrix ( List<scalar> list )
	{
	    Eigen::MatrixXd matrix(list.size(),1);
	    for(label i = 0; i < matrix.rows(); i++)
	    {
	        matrix(i,0) = list[i];
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
                scalar faceY =
                    mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].y();
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
                    gTrue[timeI][faceI] = heatFlux_steady[faceI] * (1 + 0.5 * (std::sin ( 2 * M_PI * (maxFrequency / endTime * timeSteps[timeI]) * timeSteps[timeI])));
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
                    //ITHACAutilities::assignBC(T, patchI, - Tf / k);
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
	    Info << "Thermal diffusivity = " << diffusivity << " m2/s"<< endl;

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
	    Ttrue.append(T);
            while (runTime.loop())
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;
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
                                             std::to_string(timeSteps[timeI]), directOutputFolder,
                                             "gTrue");
	        Ttrue.append(T);
	        Info << endl;
                
                runTime.printExecutionTime(Info);
                runTime.write();
            }
	    
	    Tmeas = fieldValueAtThermocouples(Ttrue);
	    forAll(samplingTime, sampleI)
	    {
	        std::cout << "debug: Tmeas = " << Tmeas(sampleI * thermocouplesNum + 5) << std::endl;
	    }

            Info<< "End\n" << endl;
	}

        void inverseProblemPostProcess(word outputFolder)
        {
            fvMesh& mesh(_mesh());
            Eigen::VectorXd heatFluxL2norm(timeSteps.size());
            Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());
            Info << "debug : g.size() = " << g.size() << endl;
            Info << "debug : timeSteps.size() = " << timeSteps.size() << endl;

            Foam::vector probe1 (0.0, 0.0, 0.6);
            label cellProbe1 = mesh.findCell(probe1);
            Eigen::VectorXd probe1_gTrue(timeSteps.size());
            Eigen::VectorXd probe1_gRec = probe1_gTrue;
            Eigen::VectorXd probe1_Ttrue = probe1_gTrue;
            Eigen::VectorXd probe1_Trec = probe1_gTrue;

            Foam::vector probe2 (0.0, 0.015, 0.6);
            label cellProbe2 = mesh.findCell(probe2);
            Eigen::VectorXd probe2_Ttrue(timeSteps.size());
            Eigen::VectorXd probe2_Trec = probe2_Ttrue;

            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "Treconstructed",
                                  outputFolder);
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
                probe1_gTrue(timeI) = gTrueField[cellProbe1];
                probe1_gRec(timeI) = gParametrizedField[cellProbe1];
                probe1_Ttrue(timeI) = Ttrue[timeI][cellProbe1];
                probe1_Trec(timeI) = Trec[timeI][cellProbe1];

                probe2_Ttrue(timeI) = Ttrue[timeI][cellProbe2];
                probe2_Trec(timeI) = Trec[timeI][cellProbe2];

                scalar EPS = 1e-6;
                volScalarField relativeErrorField(gDiffField);
                for (label i = 0; i < relativeErrorField.internalField().size(); i++)
                {
                    if (std::abs(gTrueField.ref()[i]) < EPS)
                    {
                        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
                    }
                    else
                    {
                        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / gTrueField.ref()[i];
                    }
                }
                ITHACAstream::exportSolution(relativeErrorField,
                                             std::to_string(timeSteps[timeI]), outputFolder,
                                             "relativeErrorField");
                heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, relativeErrorField,
                                            "hotSide");
                heatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh, relativeErrorField,
                                            "hotSide");
            }
            ITHACAstream::exportMatrix(heatFluxL2norm, "relErrL2norm", "eigen", outputFolder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "relErrLinfNorm", "eigen", outputFolder);
            Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
            Eigen::VectorXd samplingTime_eigen = Foam2Eigen::List2EigenMatrix(samplingTime);
            ITHACAstream::exportMatrix(timeSteps_eigen, "timeSteps", "eigen", outputFolder);
            ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe1_gTrue, "probe1_gTrue", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe1_gRec, "probe1_gRec", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe1_Ttrue, "probe1_Ttrue", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe1_Trec, "probe1_Trec", "eigen", outputFolder);

            ITHACAstream::exportMatrix(probe2_Ttrue, "probe2_Ttrue", "eigen", outputFolder);
            ITHACAstream::exportMatrix(probe2_Trec, "probe2_Trec", "eigen", outputFolder);
        }


	//void inverseProblemPostProcess(word outputFolder)
	//{
	//    fvMesh& mesh(_mesh());
	//    Eigen::VectorXd heatFluxL2norm(timeSteps.size());
	//    Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());
	//    forAll(timeSteps, timeI)
	//    {
	//        volScalarField gTrueField = list2Field(gTrue[timeI]);
	//	volScalarField gParametrizedField = list2Field(g[timeI]);
        //        volScalarField gDiffField = gParametrizedField - gTrueField;
	//	ITHACAstream::exportSolution(gTrueField,
        //                                     std::to_string(timeSteps[timeI]), outputFolder,
        //                                     "gTrue");
	//	ITHACAstream::exportSolution(Ttrue[timeI],
        //                                     std::to_string(timeSteps[timeI]), outputFolder,
        //                                     "Ttrue");
        //        ITHACAstream::exportSolution(gDiffField,
        //                                     std::to_string(timeSteps[timeI]), outputFolder,
        //                                     "gDiffField");
        //        scalar EPS = 1e-6;
        //        volScalarField relativeErrorField(gDiffField);
        //        for (label i = 0; i < relativeErrorField.internalField().size(); i++)
        //        {
        //            if (std::abs(gTrueField.ref()[i]) < EPS)
        //            {
        //                relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
        //            }
        //            else
        //            {
        //                relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / gTrueField.ref()[i];
        //            }
        //        }
        //        ITHACAstream::exportSolution(relativeErrorField,
        //                                     std::to_string(timeSteps[timeI]), outputFolder,
        //                                     "relativeErrorField");
        //        heatFluxL2norm(timeI) = ITHACAutilities::L2normOnPatch(mesh, relativeErrorField,
        //                                    "hotSide");
        //        heatFluxLinfNorm(timeI) = ITHACAutilities::LinfNormOnPatch(mesh, relativeErrorField,
        //                                    "hotSide");
	//    }
	//    ITHACAstream::exportMatrix(heatFluxL2norm, "relErrL2norm", "eigen", outputFolder);
	//    ITHACAstream::exportMatrix(heatFluxLinfNorm, "relErrLinfNorm", "eigen", outputFolder);
	//    Eigen::VectorXd timeSteps_eigen = Foam2Eigen::List2EigenMatrix(timeSteps);
	//    Eigen::VectorXd samplingTime_eigen = Foam2Eigen::List2EigenMatrix(samplingTime);
	//    ITHACAstream::exportMatrix(timeSteps_eigen, "timeSteps", "eigen", outputFolder);
	//    ITHACAstream::exportMatrix(samplingTime_eigen, "samplingTime", "eigen", outputFolder);
	//}
};

#endif

