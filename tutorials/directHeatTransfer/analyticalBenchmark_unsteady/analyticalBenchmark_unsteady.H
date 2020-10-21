class analyticalBenchmark_unsteady: public laplacianProblem
{
    public:
        analyticalBenchmark_unsteady(int argc, char* argv[])
            :
            DT("DT", dimensionSet(1, 1, -3, -1, 0, 0, 0), 1.0)
        {
            _args = autoPtr<argList>
                    (
                        new argList(argc, argv)
                    );

            if (!_args->checkRootCase())
            {
                Foam::FatalError.exit();
            }

            argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
            _simple = autoPtr<simpleControl>
                      (
                          new simpleControl
                          (
                              mesh
                          )
                      );
#include "createFields.H"
#include "createFvOptions.H"
            ITHACAdict = new IOdictionary
            (
                IOobject
                (
                    "ITHACAdict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            para = ITHACAparameters::getInstance(mesh, runTime);
            offline = ITHACAutilities::check_off();
            podex = ITHACAutilities::check_pod();
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
        }

        ITHACAparameters* para;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// Dummy thermal conductivity with unitary value
        dimensionedScalar DT;

        double k;  //thermal diffusivity [W/(m K)]
        double H;  //heat transfer coefficient [W/(m2 K)]

        //Boudary Conditions
        List<scalar> heatFlux_gammaEx1;
        List<scalar> heatFlux_gammaEx2;
        List<scalar> heatFlux_gammaEx3;
        List<scalar> heatFlux_gammaEx4;
        List<scalar> heatFlux_hotSide;
        label hotSide_ind;
        label coldSide_ind;
        label gammaEx1_ind;
        label gammaEx2_ind;
        label gammaEx3_ind;
        label gammaEx4_ind;

        PtrList<volScalarField> gField;
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;

        List<scalar> Tf;
        List<scalar> refGrad;
        List<scalar> valueFraction;
        List<scalar> homogeneousBCcoldSide;

        double a = 5;
        double b = 10;
        double c = 15;
        double d = 20;

        scalar rho;
        scalar Cp;

        void set_valueFraction()
        {
            fvMesh& mesh = _mesh();
            valueFraction.resize(mesh.boundaryMesh()["coldSide"].size());
            homogeneousBCcoldSide.resize(mesh.boundaryMesh()["coldSide"].size());
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            Eigen::VectorXd faceCellDist =
                ITHACAutilities::boudaryFaceToCellDistance(mesh, coldSide_ind);
            forAll (valueFraction, faceI)
            {
                scalar faceDist = faceCellDist(faceI);
                valueFraction[faceI] =  1.0 / (1.0 + (k / H / faceDist));
                homogeneousBCcoldSide[faceI] =  0;
            }
            refGrad = homogeneousBCcoldSide;
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
                Tf[faceI] =  k / H * b * faceY + (a + b + c) * diffusivity.value() * time + a /
                             2 * faceX * faceX +
                             b / 2 * faceY * faceY + c / 2 * faceZ * faceZ + d;

                if (faceI == 0 )
                {
                    Info << "debug :" << endl;
                    Info << Tf[faceI] << "\n" << endl;
                }
            }
            List<scalar> RobinBC = Tf;
            ITHACAutilities::assignMixedBC(T, mesh.boundaryMesh().findPatchID("coldSide"),
                                           RobinBC, refGrad, valueFraction);
        }

        void solveUnsteady()
        {
            //restart();
            volScalarField& T(_T());
            volScalarField T_true(T);
            Foam::Time& runTime = _runTime();
            word outputFolder = "./ITHACAoutput/test/";
            dimensionedScalar diffusivity("diffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0),
                                          k / (rho * Cp));
            Info << "Thermal diffusivity = " << diffusivity << " m2/s" << endl;
            simpleControl& simple = _simple();
            fv::options& fvOptions(_fvOptions());
            assignBC();
            assignIF();
            List<scalar> temperatureOverTime;
            Eigen::VectorXd relErr_L2;
            Eigen::VectorXd relErr_Linf;
            Eigen::VectorXd timeVec;
            int timeI = 0;
            scalar relErrL2_time = 0;
            scalar relErrLinf_time = 0;

            while (runTime.loop())
            {
                scalar time = runTime.value();
                Info << "Time = " << runTime.timeName() << nl << endl;
                updateBC(diffusivity);

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
                                             outputFolder,
                                             "T");

                // setting analytical solution
                for (label i = 0; i < T_true.internalField().size(); i++)
                {
                    auto cx = T_true.mesh().C()[i].component(vector::X);
                    auto cy = T_true.mesh().C()[i].component(vector::Y);
                    auto cz = T_true.mesh().C()[i].component(vector::Z);
                    T_true.ref()[i] = (a + b + c) * diffusivity.value() * time +
                                      (a / 2 * cx * cx + b / 2 * cy * cy + c / 2 * cz * cz) + d;
                }

                Info << "Exporting analytical solution" << endl;
                scalar EPS = 1e-6;
                volScalarField TdiffField = T - T_true;
                volScalarField relativeErrorField(T_true);

                for (label i = 0; i < relativeErrorField.internalField().size(); i++)
                {
                    if (std::abs(T_true.ref()[i]) < EPS)
                    {
                        relativeErrorField.ref()[i] = (std::abs(TdiffField.ref()[i])) / EPS;
                    }
                    else
                    {
                        relativeErrorField.ref()[i] = (std::abs(TdiffField.ref()[i])) / T_true.ref()[i];
                    }
                }

                ITHACAstream::exportSolution(T_true, runTime.timeName(),
                                             outputFolder,
                                             "Ttrue");
                volScalarField error = analyticalErrorField(T, T_true);
                ITHACAstream::exportSolution(error, runTime.timeName(), outputFolder, "error");
                relErr_L2.conservativeResize(timeI + 1);
                relErr_Linf.conservativeResize(timeI + 1);
                relErr_L2(timeI) = ITHACAutilities::L2Norm(
                                       relativeErrorField);
                relErr_Linf(timeI) = ITHACAutilities::LinfNorm(relativeErrorField);
                timeVec.conservativeResize(timeI + 1);
                relErrL2_time += relErr_L2(timeI) * runTime.deltaTValue();
                if(relErr_Linf(timeI) > relErrLinf_time)
                {
                    relErrLinf_time = relErr_Linf(timeI);
                }

                timeVec(timeI) = runTime.value();
                Info << "L2 norm of the relative error = " << relErr_L2(timeI) << endl;
                Info << "Linf norm of the relative error = " <<
                     relErr_Linf(timeI) << endl;
                Info << "L2 norm of the error = " << ITHACAutilities::L2Norm(error) << endl;
                Info << "Linf norm of the error = " << ITHACAutilities::LinfNorm(error) << endl;
                runTime.printExecutionTime(Info);
                runTime.write();
                timeI++;
            }
            Info << "\nRelative error L2 norm in time = " << relErrL2_time << endl; 
            Info << "Relative error Linf norm in time = " << relErrLinf_time << endl << endl; 
            std::ofstream relErrLinf_file ("relErrLinf_time.txt");
            relErrLinf_file << std::fixed << std::setprecision(8) << relErrLinf_time;
            relErrLinf_file.close();
            std::ofstream relErrL2_file ("relErrL2_time.txt");
            relErrL2_file << std::fixed << std::setprecision(8) << relErrL2_time;
            relErrL2_file.close();


            ITHACAstream::exportMatrix(relErr_L2, "relErr_L2", "eigen", outputFolder);
            ITHACAstream::exportMatrix(relErr_Linf, "relErr_Linf", "eigen", outputFolder);
            ITHACAstream::exportMatrix(timeVec, "timeVec", "eigen", outputFolder);
            Info << "End\n" << endl;
        }
};
