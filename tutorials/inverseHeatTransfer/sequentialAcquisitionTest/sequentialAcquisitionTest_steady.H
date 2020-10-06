#ifndef sequentialAcquisitionTest_steady_H
#define sequentialAcquisitionTest_steady_H


class sequentialAcquisitionTest_steady: public inverseLaplacianProblem
{
    public:
        explicit sequentialAcquisitionTest_steady(int argc, char* argv[])
            :
            inverseLaplacianProblem(argc, argv),
            T(_T()),
            lambda(_lambda()),
            deltaT(_deltaT()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            interpolationPlaneDefined = 0;
            cgIter = 0;
            thermocouplesRead = 0;
	    folderOffline = "./ITHACAoutput/offlineParamBC_steady/";
        }
        volScalarField& T;
        volScalarField& lambda;
        volScalarField& deltaT;
        fvMesh& mesh;
        Time& runTime;

        double a;
        double b;
        double c;
        double d;

        //Boudary Conditions
        List<scalar> heatFlux_hotSide;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TredField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> lambdaFullField;
        PtrList<volScalarField> deltaTfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;


        volScalarField analyticalErrorField (volScalarField& field1,
                                             volScalarField& field2)
        {
            volScalarField fieldDiff(field1 - field2);
            return fieldDiff;
        }

        void computeTdiff()
        {
            TdiffField = TredField;
            forAll(TdiffField, testI)
            {
                TdiffField[testI] = TredField[testI] - TfullField[0];
            }
        }

        int isInPlane(double cx, double cz)
        {
            return (cx >= interpolationPlane.minX -
                    interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz >= interpolationPlane.minZ - interpolationPlane.thermocoupleCellDim[2] / 4 &&
                    cx <= interpolationPlane.maxX + interpolationPlane.thermocoupleCellDim[0] / 4 &&
                    cz <= interpolationPlane.maxZ + interpolationPlane.thermocoupleCellDim[2] / 4
                   );
        }

        void postProcess(word folder, List<scalar> samplingTime, word heatFluxFieldName, scalar innerField = 0.0)
        {
            Info << "Computing errors" << endl;
            gFullField.resize(0);
            ITHACAstream::read_fields(gFullField, heatFluxFieldName,
                                      folder);
            PtrList<volScalarField> gTrueField;
            ITHACAstream::read_fields(gTrueField, "gTrue",
                                      folder);
            label Nsolutions = gFullField.size();
            Eigen::MatrixXd heatFluxL2norm;
            heatFluxL2norm.resize(Nsolutions, 1);
            Eigen::MatrixXd heatFluxLinfNorm = heatFluxL2norm;
            forAll(samplingTime, solutionI)
            {
                volScalarField gDiffField = gFullField[solutionI] - gTrueField[solutionI];
                ITHACAstream::exportSolution(gDiffField,
                                             std::to_string(samplingTime[solutionI]), folder,
                                             "gDiffField");

	        scalar EPS = 1e-6;
	        volScalarField relativeErrorField(gDiffField);
	        for (label i = 0; i < relativeErrorField.internalField().size(); i++)
                {
		    if (std::abs(gTrueField[solutionI].ref()[i]) < EPS)
		    {
		        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / EPS;
		    }
	            else
	            {
                        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) / gTrueField[solutionI].ref()[i];
	            }
                }
                ITHACAstream::exportSolution(relativeErrorField,
                                             std::to_string(samplingTime[solutionI]), folder,
                                             "relativeErrorField");
                heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh, relativeErrorField,
                                            "hotSide");
                heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh, relativeErrorField,
                                            "hotSide");

            }
            ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                                       folder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                                       folder);
        }

	void set_Tf(List<scalar> _Tf)
        {
            fvMesh& mesh = _mesh();
	    volScalarField& T(_T());
            // Setting BC at the cold side
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            label coldSideSize = T.boundaryField()[coldSide_ind].size();
            Tf = _Tf;
            refGrad.resize(coldSideSize);
            valueFraction.resize(coldSideSize);
	}

        void solveTrue()
        {
            volScalarField& T = _T();
            fvMesh& mesh = _mesh();
            Foam::Time& runTime = _runTime();
            //ITHACAutilities::assignIF(T, homogeneousBC);
            simpleControl simple(mesh);
            Info << "DT = " << DT << endl;
	    assignDirectBC();

            while (simple.loop())
            {
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        fvm::laplacian(DT, T)
                    );
                    TEqn.solve();
                }
            }
            ITHACAstream::exportSolution(T,
                                         "1", "./ITHACAoutput/debug",
                                         "T0");
	}

};

#endif
