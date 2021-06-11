class unsteadyTest: public inverseLaplacianProblem
{
    public:
        explicit unsteadyTest(int argc, char* argv[])
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
	
	scalar rho;
	scalar Cp;
	scalar alpha;

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
            }
	}

	void updateBC(dimensionedScalar diffusivity)
	{
	    Info << "Updating BC" << endl;
            Foam::Time& runTime = _runTime();
	    scalar time = runTime.value();
            fvMesh& mesh = _mesh();
            volScalarField& T(_T());
	}

	void solveUnsteady()
	{
            //restart();
            volScalarField& T(_T());
            volScalarField T_true(T);
            Foam::Time& runTime = _runTime();
	    word folderOffline = "./ITHACAoutput/test/";
	    dimensionedScalar diffusivity("diffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0), alpha);
	    Info << "Thermal diffusivity = " << diffusivity << " m2/s"<< endl;

            simpleControl& simple = _simple();
	    fv::options& fvOptions(_fvOptions());

	    assignBC();

	    List<scalar> temperatureOverTime;

            while (runTime.loop())
            {
		scalar time = runTime.value();
                Info<< "Time = " << runTime.timeName() << nl << endl;

                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::ddt(T) - fvm::laplacian(diffusivity, T)
                    );

                    TEqn.solve();
                }
                ITHACAstream::exportSolution(T, runTime.timeName(),
                                 folderOffline,
                                 "T");

                runTime.printExecutionTime(Info);
		runTime.write();
            }
            Info<< "End\n" << endl;
	}
};
