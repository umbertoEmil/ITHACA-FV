class analyticalBenchmark_unsteady: public inverseLaplacianProblem
{
    public:
        explicit analyticalBenchmark_unsteady(int argc, char* argv[])
            :
            inverseLaplacianProblem(argc, argv)
        {
	    fvMesh& mesh = _mesh();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
        }
        
	//Boudary Conditions
        List<scalar> heatFlux_gammaEx1;
        List<scalar> heatFlux_gammaEx2;
        List<scalar> heatFlux_gammaEx3;
        List<scalar> heatFlux_gammaEx4;
        List<scalar> heatFlux_hotSide;
        label gammaEx1_ind;
        label gammaEx2_ind;
        label gammaEx3_ind;
        label gammaEx4_ind;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;

	double a = 5;
        double b = 10;
        double c = 15;
        double d = 20;
	
	scalar rho;
	scalar Cp;

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

        void readThermocouples()
        {
	    Info << "WARNING: readThermocouples not implemented!" << endl;
        }

	void assignIF()
	{
            volScalarField& T(_T());
            for (label i = 0; i < T.internalField().size(); i++)
            {
                auto cx = T.mesh().C()[i].component(vector::X);
                auto cy = T.mesh().C()[i].component(vector::Y);
                auto cz = T.mesh().C()[i].component(vector::Z);
                T.ref()[i] =  (a / 2 * cx * cx + b / 2 * cy * cy + c / 2 * cz * cz) + d;
            }
	}

	void assignBC()
	{
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());

            List<scalar> RobinBC = Tf;
            set_valueFraction();
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(T, patchI, RobinBC, refGrad, valueFraction);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(T, patchI, heatFlux_hotSide);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
                {
                    ITHACAutilities::assignBC(T, patchI, heatFlux_gammaEx1);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
                {
                    ITHACAutilities::assignBC(T, patchI, heatFlux_gammaEx2);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
                {
                    ITHACAutilities::assignBC(T, patchI, heatFlux_gammaEx3);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
                {
                    ITHACAutilities::assignBC(T, patchI, heatFlux_gammaEx4);
                }
            }
	}

	void updateBC(dimensionedScalar diffusivity)
	{
	    Info << "Updating BC" << endl;
            Foam::Time& runTime = _runTime();
	    scalar time = runTime.value();
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());

            forAll(Tf, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].y();
                scalar faceZ =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();

                Tf[faceI] =  k / H * b * faceY + (a + b + c) * diffusivity.value() * time + a / 2 * faceX * faceX +
                                    b / 2 * faceY * faceY + c / 2 * faceZ * faceZ + d;
                if (faceI == 0 )
                {
                    Info << "debug :" << endl;
                    Info << Tf[faceI] << "\n"<< endl;
                }

            }
            List<scalar> RobinBC = Tf;
	    ITHACAutilities::assignMixedBC(T, mesh.boundaryMesh().findPatchID("coldSide"), RobinBC, refGrad, valueFraction);
	}

	void solveUnsteady()
	{
            //restart();
            volScalarField& T(_T());
            volScalarField T_true(T);
            Foam::Time& runTime = _runTime();
	    word folderOffline = "./ITHACAoutput/test/";
	    dimensionedScalar diffusivity("diffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0), k / (rho * Cp));
	    Info << "Thermal diffusivity = " << diffusivity << " m2/s"<< endl;

            simpleControl& simple = _simple();
	    fv::options& fvOptions(_fvOptions());

	    assignBC();
	    assignIF();

	    List<scalar> temperatureOverTime;
	    bool updatedBC = 0;

            while (runTime.loop())
            {
		scalar time = runTime.value();
                Info<< "Time = " << runTime.timeName() << nl << endl;
		if(!updatedBC)
		{
	            //updateBC(diffusivity);
		    updatedBC = 1;
		}

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
                ITHACAstream::exportSolution(T, runTime.timeName(),
                                 folderOffline,
                                 "T");

                // setting analytical solution
                for (label i = 0; i < T_true.internalField().size(); i++)
                {
                    auto cx = T_true.mesh().C()[i].component(vector::X);
                    auto cy = T_true.mesh().C()[i].component(vector::Y);
                    auto cz = T_true.mesh().C()[i].component(vector::Z);
                    T_true.ref()[i] = (a+b+c) * diffusivity.value() * time + (a / 2 * cx * cx + b / 2 * cy * cy + c / 2 * cz * cz) + d;
                }

	        Info << "Exporting analytical solution" << endl;
                ITHACAstream::exportSolution(T_true, runTime.timeName(),
                                 folderOffline,
                                 "Ttrue");
                volScalarField error = analyticalErrorField(T, T_true);
                ITHACAstream::exportSolution(error, runTime.timeName(), folderOffline, "error");
                Info << "L2 norm of the relative error = " << ITHACAutilities::error_fields(
                         T, T_true) << endl;
                Info << "Linf norm of the relative error = " <<
                     ITHACAutilities::relError_fields_LinfNorm(T, T_true) << endl;
                Info << "L2 norm of the error = " << ITHACAutilities::L2norm(error) << endl;
                Info << "Linf norm of the error = " << ITHACAutilities::LinfNorm(error) << endl;

                runTime.printExecutionTime(Info);
		runTime.write();
            }
            Info<< "End\n" << endl;
	}
};
