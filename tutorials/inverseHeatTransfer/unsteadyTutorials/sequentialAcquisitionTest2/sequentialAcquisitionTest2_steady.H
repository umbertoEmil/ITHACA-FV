#ifndef sequentialAcquisitionTest2_steady_H
#define sequentialAcquisitionTest2_steady_H

#include "inverseLaplacianProblem_paramBC.H"


class sequentialAcquisitionTest2_steady: public inverseLaplacianProblem_paramBC
{
    public:
        explicit sequentialAcquisitionTest2_steady(int argc, char* argv[])
            :
            inverseLaplacianProblem_paramBC(argc, argv),
            T(_T()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
            folderOffline = "./ITHACAoutput/offlineParamBC_steady/";
        }
        volScalarField& T;
        fvMesh& mesh;
        Time& runTime;

        double a;
        double b;
        double c;
        double d;

        //Boudary Conditions
        List<scalar> heatFlux_hotSide;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;

        void postProcess(word folder, List<scalar> samplingTime, word heatFluxFieldName,
                         Foam::vector probe1, Foam::vector probe2, scalar innerField = 0.0)
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
            label cellProbe1 = mesh.findCell(probe1);
            Eigen::VectorXd probe1_gRec(Nsolutions);
            Eigen::VectorXd probe1_Trec = probe1_gRec;
            label cellProbe2 = mesh.findCell(probe2);
            Eigen::VectorXd probe2_Trec(Nsolutions);
            PtrList<volScalarField> Trec;
            ITHACAstream::read_fields(Trec, "T",
                                      folder);
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
                        relativeErrorField.ref()[i] = (std::abs(gDiffField.ref()[i])) /
                                                      gTrueField[solutionI].ref()[i];
                    }
                }

                ITHACAstream::exportSolution(relativeErrorField,
                                             std::to_string(samplingTime[solutionI]), folder,
                                             "relativeErrorField");
                heatFluxL2norm(solutionI) = ITHACAutilities::L2normOnPatch(mesh,
                                            relativeErrorField,
                                            "hotSide");
                heatFluxLinfNorm(solutionI) = ITHACAutilities::LinfNormOnPatch(mesh,
                                              relativeErrorField,
                                              "hotSide");
                probe1_gRec(solutionI) = gFullField[solutionI][cellProbe1];
                probe1_Trec(solutionI) = Trec[solutionI][cellProbe1];
                probe2_Trec(solutionI) = Trec[solutionI][cellProbe2];
            }
            ITHACAstream::exportMatrix(heatFluxL2norm, "relError_L2norm", "eigen",
                                       folder);
            ITHACAstream::exportMatrix(heatFluxLinfNorm, "relError_LinfNorm", "eigen",
                                       folder);
            ITHACAstream::exportMatrix(probe1_gRec, "probe1_gRec", "eigen", folder);
            ITHACAstream::exportMatrix(probe1_Trec, "probe1_Trec", "eigen", folder);
            ITHACAstream::exportMatrix(probe2_Trec, "probe2_Trec", "eigen", folder);
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
