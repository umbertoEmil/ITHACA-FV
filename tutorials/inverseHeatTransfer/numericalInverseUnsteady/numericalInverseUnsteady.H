class numericalInverseUnsteady: public inverseHeatTransferProblem
{
    public:
        explicit numericalInverseUnsteady(int argc, char* argv[], scalar _diffusivity)
            :
            inverseHeatTransferProblem(argc, argv, _diffusivity)
        {
	    fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;

	    T0 = 0;//[K]
            
        }
	word directOutputFolder = "./ITHACAoutput/directProblem/";
	scalar a = 0;
	scalar b = 0;
	scalar c = 0;
	scalar d = 0;
	PtrList<volScalarField> Ttrue;

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
                auto cx = T.mesh().C()[i].component(vector::X);
                auto cy = T.mesh().C()[i].component(vector::Y);
                auto cz = T.mesh().C()[i].component(vector::Z);
		T.ref()[i] =  T0;
            }
        }
	void assignT0_IF(volScalarField& T0_field)
	{
	    assignTrueIF();
            volScalarField& T(_T());
	    ITHACAutilities::assignIF(T0_field, T);
	}

        void set_gTrue()
        {
	    fvMesh& mesh(_mesh());
            volScalarField& T(_T());
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
            gTrue.resize(timeSteps.size());

	    forAll(timeSteps, timeI)
	    {
		// Setting gTrue 
                gTrue[timeI].resize(hotSideSize);
                forAll(gTrue[timeI], faceI)
                {
                    scalar faceX =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].x();
                    scalar faceY =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].y();
                    scalar faceZ =
                        mesh.boundaryMesh()[hotSide_ind].faceCentres()[faceI].z();
                    gTrue[timeI][faceI] = a * faceX + b * faceZ + c * timeSteps[timeI] + d;
                }
	    }
        }

        void assignDirectBC(label timeI)
	{
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());
            // Setting BC at the cold side
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            label coldSideSize = T.boundaryField()[coldSide_ind].size();
            Tf.resize(coldSideSize);
            forAll(Tf, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();

                Tf[faceI] = (350 - 30 * faceZ);
		Tf[faceI] -= Tf[faceI] * 0.1 * mesh.time().value();
            }
	    
	    set_valueFraction();
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
                        fvm::ddt(T) - fvm::laplacian(diffusivity, T)
                    );

                    fvOptions.constrain(TEqn);
                    TEqn.solve();
                    fvOptions.correct(T);
                }
                ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI]),
                                 directOutputFolder,
                                 "Ttrue");
	        Ttrue.append(T);
		Info << endl;
                
                runTime.printExecutionTime(Info);
                runTime.write();
            }
	    Tmeas = fieldValueAtThermocouples(Ttrue);

            Info<< "End\n" << endl;
	}

	void inverseProblemPostProcess(word outputFolder)
	{
	    fvMesh& mesh(_mesh());
	    Eigen::VectorXd heatFluxL2norm(timeSteps.size());
	    Eigen::VectorXd heatFluxLinfNorm(timeSteps.size());
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
                scalar EPS = 1e-6;
                volScalarField relativeErrorField(gDiffField);
                for (label i = 0; i < relativeErrorField.internalField().size(); i++)
                {
                    if (gTrueField.ref()[i] < EPS)
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
	}
};
