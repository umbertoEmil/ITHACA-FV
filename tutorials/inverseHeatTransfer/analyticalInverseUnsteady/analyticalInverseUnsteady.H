class analyticalInverseUnsteady: public inverseHeatTransferProblem
{
    public:
        explicit analyticalInverseUnsteady(int argc, char* argv[], scalar _diffusivity)
            :
            inverseHeatTransferProblem(argc, argv, _diffusivity)
        {
	    fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;

	    T0 = 30;//[K]
            
        }
	word directOutputFolder = "./ITHACAoutput/directProblem/";
	scalar a = 5;
	scalar b = 1;
	scalar c = 15;
	scalar d = 20;
	unsigned exponentialSolution;
	unsigned quadraticSolution;
	PtrList<volScalarField> Tanalytical;
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

        volScalarField analyticalErrorField (volScalarField& field1,
                                             volScalarField& field2)
        {
            volScalarField fieldDiff(field1 - field2);
            return fieldDiff;
        }

        void assignTrueIF()
        {
            volScalarField& T(_T());
            for (label i = 0; i < T.internalField().size(); i++)
            {
                auto cx = T.mesh().C()[i].component(vector::X);
                auto cy = T.mesh().C()[i].component(vector::Y);
                auto cz = T.mesh().C()[i].component(vector::Z);
		if(quadraticSolution)
		{
		    T.ref()[i] =  a * cx * cx + b * cy * cy + c * cz * cz;
		}
		else if(exponentialSolution)
		{
                    T.ref()[i] =  a * cx * cx + b * cx * cz * std::exp(cy / std::sqrt(diffusivity.value())) + c * cz *cz + d;
		}
            }
        }
	void assignT0_IF(volScalarField& T0_field)
	{
	    assignTrueIF();
            volScalarField& T(_T());
	    ITHACAutilities::assignIF(T0_field, T);
	}

        void analyticalSolution()
        {
	    fvMesh& mesh(_mesh());
            volScalarField& T(_T());
	    Tanalytical.resize(0);
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            label hotSideSize = T.boundaryField()[hotSide_ind].size();
            gTrue.resize(timeSteps.size());

	    forAll(timeSteps, timeI)
	    {
                for (label i = 0; i < T.internalField().size(); i++)
                {
                    auto cx = T.mesh().C()[i].component(vector::X);
                    auto cy = T.mesh().C()[i].component(vector::Y);
                    auto cz = T.mesh().C()[i].component(vector::Z);
		    if(quadraticSolution)
		    {
                        T.ref()[i] =  a * cx * cx + b * cy * cy + c * cz *cz 
			    + diffusivity.value() * 2 * (a + b + c) * timeSteps[timeI];
		    }
		    else if(exponentialSolution)
		    {
                        T.ref()[i] =  a * cx * cx + b * cx * cz * std::exp(timeSteps[timeI] + 
		            cy / std::sqrt(diffusivity.value())) + c * cz * cz + d +
	                    diffusivity.value() * 2 * (a + c) * timeSteps[timeI];
		    }
                }
		Tanalytical.append(T);
                ITHACAstream::exportSolution(T, std::to_string(timeSteps[timeI]),
                                 directOutputFolder,
                                 "Tanalytical");
	        
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
		    if(quadraticSolution)
		    {
                        gTrue[timeI][faceI] = k * 2 * b * faceY;
		    }
		    else if(exponentialSolution)
		    {
		        gTrue[timeI][faceI] = k * b / std::sqrt(diffusivity.value()) * faceX * faceZ *
	                    std::exp(timeSteps[timeI] + faceY / std::sqrt(diffusivity.value()));
		    }
                }
	    }

	    Tmeas = fieldValueAtThermocouples(Tanalytical);
            std::cout << "debug: Tmeas true = " << Tmeas << std::endl;
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
		if(quadraticSolution)
		{
		    //Tf[faceI] = - k * 2 * b * faceY;
                    Tf[faceI] = (k / H) * 2 * b * faceY + (a * faceX * faceX +
		        b * faceY + c * faceZ * faceZ
			+ 2 * diffusivity.value() * ( a + b + c ) * mesh.time().value());
		}
		else if(exponentialSolution)
		{
		    //Tf[faceI] = - k * b / std::sqrt(diffusivity.value()) *
                    //    (faceX * faceZ * std::exp(mesh.time().value() +
                    //    faceY / std::sqrt(diffusivity.value())));
		    Tf[faceI] = (k * b) / (H * std::sqrt(diffusivity.value())) *  
		        (faceX * faceZ * std::exp(mesh.time().value() + 
		        faceY / std::sqrt(diffusivity.value()))) + (a * faceX * faceX +
		        b * faceX * faceZ * std::exp(faceY / std::sqrt(diffusivity.value()) 
		        + mesh.time().value()) + c * faceZ * faceZ + d +
		        2 * diffusivity.value() * (a + c) * mesh.time().value());
		}
            }

            // Setting BC at gammaEx1 
            label gammaEx1_ind = mesh.boundaryMesh().findPatchID("gammaEx1");
            label gammaEx1Size = T.boundaryField()[gammaEx1_ind].size();
	    List<scalar> heatFlux_gammaEx1;
            heatFlux_gammaEx1.resize(gammaEx1Size);
            forAll(heatFlux_gammaEx1, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx1_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx1_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx1_ind].faceCentres()[faceI].z();
                
		if(quadraticSolution)
		{
                    heatFlux_gammaEx1[faceI] = - k * 2 * c * faceZ;
		}
		else if(exponentialSolution)
		{
                    heatFlux_gammaEx1[faceI] = - k * (b * faceX * 
		        std::exp(faceY / std::sqrt(diffusivity.value()) + 
		        mesh.time().value()) + 2 * c * faceZ);
		}
            }
            
            // Setting BC at gammaEx2 
            label gammaEx2_ind = mesh.boundaryMesh().findPatchID("gammaEx2");
            label gammaEx2Size = T.boundaryField()[gammaEx2_ind].size();
	    List<scalar> heatFlux_gammaEx2;
            heatFlux_gammaEx2.resize(gammaEx2Size);
            forAll(heatFlux_gammaEx2, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].z();

                if(quadraticSolution)
                {
                    heatFlux_gammaEx2[faceI] = - k * 2 * a * faceX; 
                }
                else if(exponentialSolution)
                {
                    heatFlux_gammaEx2[faceI] = - k * (2 * a * faceX + b * faceZ * 
		        std::exp(mesh.time().value() + faceY / std::sqrt(diffusivity.value())));
                }
            }
            
            // Setting BC at gammaEx3 
            label gammaEx3_ind = mesh.boundaryMesh().findPatchID("gammaEx3");
            label gammaEx3Size = T.boundaryField()[gammaEx3_ind].size();
	    List<scalar> heatFlux_gammaEx3;
            heatFlux_gammaEx3.resize(gammaEx3Size);
            forAll(heatFlux_gammaEx3, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx3_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx3_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx3_ind].faceCentres()[faceI].z();
                if(quadraticSolution)
                {
                    heatFlux_gammaEx3[faceI] = k * 2 * c * faceZ;
                }
                else if(exponentialSolution)
                {
                    heatFlux_gammaEx3[faceI] = k * (b * faceX * 
		        std::exp(faceY / std::sqrt(diffusivity.value()) + 
		        mesh.time().value()) + 2 * c * faceZ);
                }
            }
            
            // Setting BC at gammaEx4 
            label gammaEx4_ind = mesh.boundaryMesh().findPatchID("gammaEx4");
            label gammaEx4Size = T.boundaryField()[gammaEx4_ind].size();
	    List<scalar> heatFlux_gammaEx4;
            heatFlux_gammaEx4.resize(gammaEx4Size);
            forAll(heatFlux_gammaEx4, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].z();
                if(quadraticSolution)
                {
                    heatFlux_gammaEx4[faceI] = k * (2 * a * faceX);
                }
                else if(exponentialSolution)
                {
                    heatFlux_gammaEx4[faceI] = k * (2 * a * faceX +
		        b * faceZ * std::exp(faceY / std::sqrt(diffusivity.value()) +
		        mesh.time().value()));
                }

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
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
                {   
                    ITHACAutilities::assignBC(T, patchI, - heatFlux_gammaEx1 / k);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
                {   
                    ITHACAutilities::assignBC(T, patchI, - heatFlux_gammaEx2 / k);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
                {   
                    ITHACAutilities::assignBC(T, patchI, - heatFlux_gammaEx3 / k);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
                {   
                    ITHACAutilities::assignBC(T, patchI, - heatFlux_gammaEx4 / k);
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

            Info<< "End\n" << endl;
	}

	void directProblemPostProcess()
	{
	    Eigen::VectorXd relError_L2norm;
	    relError_L2norm.resize(timeSteps.size());
	    Eigen::VectorXd relError_LinfNorm = relError_L2norm;
	    Eigen::VectorXd Error_L2norm = relError_L2norm;
	    Eigen::VectorXd Error_LinfNorm = relError_L2norm;
	    Eigen::VectorXd timeEigen = relError_L2norm;


	    forAll(timeSteps, timeI)
	    {
                volScalarField error = Tanalytical[timeI] - Ttrue[timeI];
                for (label i = 0; i < error.internalField().size(); i++)
                {
                    error.ref()[i] =  std::abs(error.ref()[i]); 
                }

                timeEigen(timeI) = timeSteps[timeI];
                ITHACAstream::exportSolution(error, std::to_string(timeSteps[timeI]),
		    directOutputFolder, "directError");
		relError_L2norm(timeI) = ITHACAutilities::error_fields(
                         Tanalytical[timeI], Ttrue[timeI]);
	        relError_LinfNorm(timeI) = ITHACAutilities::relError_fields_LinfNorm(Tanalytical[timeI], Ttrue[timeI]);
		Error_L2norm(timeI) = ITHACAutilities::L2norm(error);
		Error_LinfNorm(timeI) = ITHACAutilities::LinfNorm(error);
	    }
	    ITHACAstream::exportMatrix(timeEigen, "timeSteps", "eigen", directOutputFolder);
	    ITHACAstream::exportMatrix(relError_L2norm, "relError_L2norm", "eigen", directOutputFolder);
	    ITHACAstream::exportMatrix(relError_LinfNorm, "relError_LinfNorm", "eigen", directOutputFolder);
	    ITHACAstream::exportMatrix(Error_L2norm, "Error_L2norm", "eigen", directOutputFolder);
	    ITHACAstream::exportMatrix(Error_LinfNorm, "Error_LinfNorm", "eigen", directOutputFolder);
	}

	void inverseProblemPostProcess(word outputFolder)
	{
	    forAll(timeSteps, timeI)
	    {
	        volScalarField gTrueField = list2Field(gTrue[timeI]);
		ITHACAstream::exportSolution(gTrueField,
                                             std::to_string(timeSteps[timeI]), outputFolder,
                                             "gTrue");
	    }
	}
};
