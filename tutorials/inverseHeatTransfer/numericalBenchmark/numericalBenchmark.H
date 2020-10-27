#ifndef numericalBenchmark_H
#define numericalBenchmark_H


class numericalBenchmark: public inverseLaplacianProblem
{
    public:
        explicit numericalBenchmark(int argc, char* argv[])
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
        }
        volScalarField& T;
        volScalarField& lambda;
        volScalarField& deltaT;
        fvMesh& mesh;
        Time& runTime;

        double g_0;
        double g_Z;
        double g_X;
        double Tf_0;
        double Tf_delta;
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

        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
            Info << "Computing errors" << endl;
            gFullField.resize(0);
            ITHACAstream::read_fields(gFullField, heatFluxFieldName,
                                      folder);
            volScalarField gTrueField = list2Field(gTrue, innerField);
            label Nsolutions = gFullField.size();
            Eigen::MatrixXd heatFluxL2norm;
            heatFluxL2norm.resize(Nsolutions, 1);
            Eigen::MatrixXd heatFluxLinfNorm = heatFluxL2norm;
            forAll(gFullField, solutionI)
            {
                volScalarField gDiffField = gFullField[solutionI] - gTrueField;
                ITHACAstream::exportSolution(gDiffField,
                                             std::to_string(solutionI + 1), folder,
                                             "gDiffField");
                ITHACAstream::exportSolution(gTrueField,
                                             std::to_string(solutionI + 1), folder,
                                             "gTrueField");

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
                                             std::to_string(solutionI + 1), folder,
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

	void set_Tf()
        {
            fvMesh& mesh = _mesh();
	    volScalarField& T(_T());
            // Setting BC at the cold side
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            label coldSideSize = T.boundaryField()[coldSide_ind].size();
            Tf.resize(coldSideSize);
            refGrad.resize(coldSideSize);
            valueFraction.resize(coldSideSize);
            forAll(Tf, faceI)
            {
                scalar faceZ =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].z();
                scalar faceY =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].y();
                scalar faceX =
                    mesh.boundaryMesh()[coldSide_ind].faceCentres()[faceI].x();

                Tf[faceI] = (Tf_delta * (1.2 - faceZ) + Tf_0);
            }
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
            word outputFolder = "./ITHACAoutput/true/";

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

            Info << "Exporting" << endl;
            ITHACAstream::exportSolution(T, "1", outputFolder, T.name());
            volScalarField gTrueField = list2Field(g);
            ITHACAstream::exportSolution(gTrueField,
                                     "1", outputFolder,
                                     "gTrue");


            //Reinitializing runTime
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 0);
        }
};

#endif
