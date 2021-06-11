#ifndef analyticalBenchmark_H
#define analyticalBenchmark_H


class analyticalBenchmark: public inverseLaplacianProblem
{
    public:
        explicit analyticalBenchmark(int argc, char* argv[])
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

        double a;
        double b;
        double c;
        double d;

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

        // To be added to ITHACAutilities
        scalar list2norm(List<scalar> list)
        {
            scalar norm = 0;
            forAll(list, I)
            {
                norm += list[I] * list[I];
            }
            return Foam::sqrt(norm);
        }

        void solveAdditional()
	{
            restart();
            fvMesh& mesh = _mesh();
            simpleControl& simple = _simple();
            volScalarField Tad(_T);
            Foam::Time& runTime = _runTime();
            set_valueFraction();
            List<scalar> RobinBC = - Tf;

            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(Tad, patchI, RobinBC, refGrad, valueFraction);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, homogeneousBC);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx1);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx2);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx3);
                }
                else if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
                {
                    ITHACAutilities::assignBC(Tad, patchI, -heatFlux_gammaEx4);
                }
            }
            #if OFVER == 6
            
                while (simple.loop(runTime))
            #else
                while (simple.loop())
            #endif
                {
                    while (simple.correctNonOrthogonal())
                    {
                        fvScalarMatrix TEqn
                        (
                            fvm::laplacian(DT, Tad)
                        );
                        TEqn.solve();
                    }
                }
		Tad_base.append(Tad);
            
                //Reinitializing runTime
                instantList Times = runTime.times();
                runTime.setTime(Times[1], 1);
                addSol = fieldValueAtThermocouples(Tad);
                ITHACAstream::exportSolution(Tad, "1",
                                             "ITHACAoutput/offlineParamBC/",
                                             "Tad");
	}


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
            Info << "Computing errors" << endl;
            gFullField.resize(0);
            ITHACAstream::read_fields(mesh, gFullField, heatFluxFieldName,
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
                //heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh, gDiffField,
                //                            "hotSide") / ITHACAutilities::L2normOnPatch(mesh, gTrueField, "hotSide");
                //heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh, gDiffField,
                //                              "hotSide") / ITHACAutilities::LinfNormOnPatch(mesh, gTrueField, "hotSide");

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

	void assignDirectBC()
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
                Tf[faceI] = k / H * (b * faceX + c) + a * faceX * faceX +
                                    b * faceX * faceY + c * faceY - a * faceZ * faceZ + c;
            }
            // Setting BC at gammaEx1
            gammaEx1_ind = mesh.boundaryMesh().findPatchID("gammaEx1");
            label gammaEx1Size = T.boundaryField()[gammaEx1_ind].size();
            heatFlux_gammaEx1.resize(gammaEx1Size);
            forAll(heatFlux_gammaEx1, faceI)
            {
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx1_ind].faceCentres()[faceI].z();
                heatFlux_gammaEx1[faceI] = k * 2 * a * faceZ;
            }
            // Setting BC at gammaEx2
            gammaEx2_ind = mesh.boundaryMesh().findPatchID("gammaEx2");
            label gammaEx2Size = T.boundaryField()[gammaEx2_ind].size();
            heatFlux_gammaEx2.resize(gammaEx2Size);
            forAll(heatFlux_gammaEx2, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx2_ind].faceCentres()[faceI].y();
                heatFlux_gammaEx2[faceI] = - k *(2 * a * faceX + b * faceY);
            }
            // Setting BC at gammaEx3
            gammaEx3_ind = mesh.boundaryMesh().findPatchID("gammaEx3");
            label gammaEx3Size = T.boundaryField()[gammaEx3_ind].size();
            heatFlux_gammaEx3.resize(gammaEx3Size);
            forAll(heatFlux_gammaEx3, faceI)
            {
                scalar faceZ =
                    mesh.boundaryMesh()[gammaEx3_ind].faceCentres()[faceI].z();
                heatFlux_gammaEx3[faceI] = 0;
            }
            // Setting BC at gammaEx4
            gammaEx4_ind = mesh.boundaryMesh().findPatchID("gammaEx4");
            label gammaEx4Size = T.boundaryField()[gammaEx4_ind].size();
            heatFlux_gammaEx4.resize(gammaEx4Size);
            forAll(heatFlux_gammaEx4, faceI)
            {
                scalar faceX =
                    mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].x();
                scalar faceY =
                    mesh.boundaryMesh()[gammaEx4_ind].faceCentres()[faceI].y();
                heatFlux_gammaEx4[faceI] = k * (b * faceY);
            }
            set_valueFraction();
            forAll(mesh.boundaryMesh(), patchI)
            {
                if (patchI == mesh.boundaryMesh().findPatchID("coldSide"))
                {
                    ITHACAutilities::assignMixedBC(T, patchI, Tf, refGrad, valueFraction);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("hotSide"))
                {
                    ITHACAutilities::assignBC(T, patchI, - g / k);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx1"))
                {
                    ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx1 / k);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx2"))
                {
                    ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx2 / k);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx3"))
                {
                    ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx3 / k);
                }

                if (patchI == mesh.boundaryMesh().findPatchID("gammaEx4"))
                {
                    ITHACAutilities::assignBC(T, patchI, -heatFlux_gammaEx4 / k);
                }
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
            ITHACAstream::exportSolution(T, "1", "./ITHACAoutput/true/", T.name());
            //Reinitializing runTime
            instantList Times = runTime.times();
            runTime.setTime(Times[1], 1);
        }


        //PtrList<volScalarField> relativeDifferenceField (PtrList<volScalarField>&
        //        field1, PtrList<volScalarField>& field2)
        //{
        //    PtrList<volScalarField> fieldDiff;
        //
        //    if (field2.size() == fieldDiff.size())
        //    {
        //        for (label i = 0; i < field1.size(); i++)
        //        {
        //            fieldDiff.append((field1[i] - field2[i]) / field2[i]);
        //        }
        //    }
        //    else if (field2.size() == 1)
        //    {
        //        for (label i = 0; i < field1.size(); i++)
        //        {
        //            fieldDiff.append((field1[i] - field2[0]) / field2[0]);
        //        }
        //    }
        //    else
        //    {
        //        Info << "WARNING:" << endl;
        //        Info << "Don't know how to behave in relativeDifferenceField" << endl;
        //        Info << "the two fields have different size" << endl;
        //        Info << "Size field1 = " << field1.size()
        //             << ", Size field2 = " << field2.size() << endl;
        //    }
        //
        //    return fieldDiff;
        //}
};

#endif
