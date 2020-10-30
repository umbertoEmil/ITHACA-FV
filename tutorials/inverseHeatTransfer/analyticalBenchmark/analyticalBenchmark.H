#ifndef analyticalBenchmark_H
#define analyticalBenchmark_H

#include "MUQ/Modeling/Distributions/Gaussian.h"

class analyticalBenchmark_paramBC: public inverseLaplacianProblem_paramBC
{
    public:
        explicit analyticalBenchmark_paramBC(int argc, char* argv[])
            :
            inverseLaplacianProblem_paramBC(argc, argv),
            T(_T()),
            mesh(_mesh()),
            runTime(_runTime())
        {
            hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
            coldSide_ind = mesh.boundaryMesh().findPatchID("coldSide");
            thermocouplesRead = 0;
        }

        volScalarField& T;
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
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;


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
                                             folderOffline,
                                             "Tad");
	}


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include"postProcess.H"
        }

	void assignDirectBC()
        {
#include"directBC.H"
	}

        void solveTrue()
        {
            word outputFolder = "./ITHACAoutput/true/";
#include"solveTrue.H"
        }

        void addNoise(double noiseStdDev_)
        {
            double variance = Tmeas.mean() * noiseStdDev_* noiseStdDev_;
            std::cout << "debug : noiseStdDev_ = " << noiseStdDev_  << std::endl;
            std::cout << "debug : Tmeas.mean() = " << Tmeas.mean()  << std::endl;
            std::cout << "debug : variance = " << variance  << std::endl;

            auto noiseDensity = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(Tmeas.size()), Eigen::MatrixXd::Identity(Tmeas.size(), Tmeas.size()) * variance);


            Tmeas += noiseDensity->Sample();
        }
};


class analyticalBenchmark_CG: public inverseLaplacianProblem_CG
{
    public:
        explicit analyticalBenchmark_CG(int argc, char* argv[])
            :
            inverseLaplacianProblem_CG(argc, argv),
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


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include"postProcess.H"
        }

	void assignDirectBC()
        {
#include"directBC.H"
	}

        void solveTrue()
        {
            word outputFolder = "./ITHACAoutput/true/";
#include"solveTrue.H"
        }

        void addNoise(double noiseStdDev_)
        {
            double variance = Tmeas.mean() * noiseStdDev_* noiseStdDev_;
            std::cout << "debug : noiseStdDev_ = " << noiseStdDev_  << std::endl;
            std::cout << "debug : Tmeas.mean() = " << Tmeas.mean()  << std::endl;
            std::cout << "debug : variance = " << variance  << std::endl;

            auto noiseDensity = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(Tmeas.size()), Eigen::MatrixXd::Identity(Tmeas.size(), Tmeas.size()) * variance);


            Tmeas += noiseDensity->Sample();

            //Set CG tolerance using Discrepancy principle
            Jtol = Tmeas.size() * variance *
                       variance / 2;
            if ( Jtol < 1e-1 )
            {
                Jtol = 1e-1;
            }
            std::cout << "debug : Tmeas2 = \n" << Tmeas.transpose() << std::endl;
        }
};
#endif
