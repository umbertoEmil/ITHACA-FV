#ifndef numericalBenchmark_H
#define numericalBenchmark_H


#include "MUQ/Modeling/Distributions/Gaussian.h"

class numericalBenchmark_CG: public inverseLaplacianProblem_CG
{
    public:
        explicit numericalBenchmark_CG(int argc, char* argv[])
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


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include "postProcess.H"
        }

	void set_Tf()
        {
#include "setTf.H"
	}

        void solveTrue()
        {
#include"solveTrue.H"
        }

        void addNoise(double noiseStdDev_)
        {
            std::cout << "debug : Tmeas1 = \n" << Tmeas.transpose() << std::endl;

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

class numericalBenchmark_paramBC: public inverseLaplacianProblem_paramBC
{
    public:
        explicit numericalBenchmark_paramBC(int argc, char* argv[])
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
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include "postProcess.H"
        }

	void set_Tf()
        {
#include "setTf.H"
	}

        void solveTrue()
        {
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

#endif
