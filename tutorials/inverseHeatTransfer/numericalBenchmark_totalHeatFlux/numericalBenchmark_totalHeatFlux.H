#ifndef numericalBenchmark_totalHeatFlux_H
#define numericalBenchmark_totalHeatFlux_H

#include "MUQ/Modeling/Distributions/Gaussian.h"

class numericalBenchmark_totalHeatFlux_paramBC: public
    inverseLaplacianProblemTotalHeatMeasure_paramBC
{
    public:
        explicit numericalBenchmark_totalHeatFlux_paramBC(int argc, char* argv[])
            :
            inverseLaplacianProblemTotalHeatMeasure_paramBC(argc, argv),
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
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;

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

        void addNoise(double noiseStdDev_, double GnoiseStdDev_)
        {
            double variance = Tmeas.mean() * noiseStdDev_ * noiseStdDev_;
            double Gvariance = std::abs(gIntegral_meas) * GnoiseStdDev_ * GnoiseStdDev_;
            std::cout << "debug : variance = " << variance  << std::endl;
            std::cout << "debug : Gvariance = " << Gvariance  << std::endl;
            auto noiseDensity = std::make_shared<muq::Modeling::Gaussian>
                                (Eigen::VectorXd::Zero(Tmeas.size()), Eigen::MatrixXd::Identity(Tmeas.size(),
                                        Tmeas.size()) * variance);
            auto GnoiseDensity = std::make_shared<muq::Modeling::Gaussian>
                                 (Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1, 1) * Gvariance);
            Tmeas += noiseDensity->Sample();
            std::cout << "debug1 G = " << gIntegral_meas << std::endl;
            gIntegral_meas += GnoiseDensity->Sample()(0, 0);
            std::cout << "debug2 G = " << gIntegral_meas << std::endl;
        }
};


class numericalBenchmark_totalHeatFlux_CG: public
    inverseLaplacianProblemTotalHeatMeasure_CG
{
    public:
        explicit numericalBenchmark_totalHeatFlux_CG(int argc, char* argv[])
            :
            inverseLaplacianProblemTotalHeatMeasure_CG(argc, argv),
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
        PtrList<volScalarField> TfullField;
        PtrList<volScalarField> TdiffField;
        PtrList<volScalarField> gFullField;
        PtrList<volScalarField> gRelErrField;


        void postProcess(word folder, word heatFluxFieldName, scalar innerField = 0.0)
        {
#include"postProcess.H"
        }

        void set_Tf()
        {
#include "setTf.H"
        }

        void solveTrue()
        {
#include"solveTrue.H"
        }


        void addNoise(double noiseStdDev_, double GnoiseStdDev_)
        {
            double variance = Tmeas.mean() * noiseStdDev_ * noiseStdDev_;
            double Gvariance = gIntegral_meas * GnoiseStdDev_ * GnoiseStdDev_;
            std::cout << "debug : variance = " << variance  << std::endl;
            std::cout << "debug : Gvariance = " << Gvariance  << std::endl;
            auto noiseDensity = std::make_shared<muq::Modeling::Gaussian>
                                (Eigen::VectorXd::Zero(Tmeas.size()), Eigen::MatrixXd::Identity(Tmeas.size(),
                                        Tmeas.size()) * variance);
            auto GnoiseDensity = std::make_shared<muq::Modeling::Gaussian>
                                 (Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1, 1) * Gvariance);
            Tmeas += noiseDensity->Sample();
            std::cout << "debug1 G = " << gIntegral_meas << std::endl;
            gIntegral_meas += GnoiseDensity->Sample()(0, 0);
            std::cout << "debug2 G = " << gIntegral_meas << std::endl;
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
