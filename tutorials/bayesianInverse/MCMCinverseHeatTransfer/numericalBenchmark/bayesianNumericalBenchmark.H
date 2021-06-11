#ifndef bayesianNumericalBenchmark_H
#define bayesianNumericalBenchmark_H


#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "numericalBenchmark.H"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/ProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"
#include "MUQ/SamplingAlgorithms/MySingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MALAProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/DRKernel.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"

class bayesianNumericalBenchmark
{
    public:
        bayesianNumericalBenchmark(numericalBenchmark& _deterministicProblem,
                                   double measurementsStdDev, double regularizationParameter, int TCrows,
                                   int TCcols, boost::property_tree::ptree pt)
            :
            deterministicProblem(_deterministicProblem),
            Nweights(deterministicProblem.gWeights.size()),
            Nthermocouples(deterministicProblem.Tmeas.size()),
            weights(Nweights),
            weightsPiece(std::make_shared<muq::Modeling::IdentityOperator>(Nweights)),
            Theta(deterministicProblem.Theta)
        {
            ThetaOp = muq::Modeling::LinearOperator::Create(Theta);
            MRFcorrelationMatrix = get_MRFcorrelationMatrix(TCrows, TCcols);
            likelihoodStdDev = measurementsStdDev;
            likelihoodCov = likelihoodStdDev * likelihoodStdDev * Eigen::MatrixXd::Identity(
                                Nthermocouples, Nthermocouples);
            likelihood = std::make_shared<muq::Modeling::Gaussian>
                         (deterministicProblem.addSol + deterministicProblem.Tmeas,
                          likelihoodCov)->AsDensity();
            priorMean = Eigen::VectorXd::Zero(Nweights);
            priorPrecision = MRFcorrelationMatrix * (2.0 * regularizationParameter /
                             (measurementsStdDev * measurementsStdDev));
            prior = std::make_shared<muq::Modeling::Gaussian>(priorMean, priorPrecision,
                    muq::Modeling::Gaussian::Mode::Precision )->AsDensity();
            prodDens = std::make_shared<muq::Modeling::DensityProduct>(2);
            graph = std::make_shared<muq::Modeling::WorkGraph>();
            assebleGraph();
            joinDens = graph->CreateModPiece("prodDens");
            problem = std::make_shared<muq::SamplingAlgorithms::SamplingProblem>(joinDens);
            mcmc = std::make_shared<muq::SamplingAlgorithms::SingleChainMCMC>(pt, problem);
        }

        numericalBenchmark& deterministicProblem;
        const unsigned int Nweights;
        const unsigned int Nthermocouples;
        Eigen::MatrixXd MRFcorrelationMatrix;

        Eigen::VectorXd weights;
        std::shared_ptr<muq::Modeling::IdentityOperator> weightsPiece;

        Eigen::MatrixXd Theta;
        std::shared_ptr<muq::Modeling::LinearOperator> ThetaOp;

        double likelihoodStdDev;
        Eigen::MatrixXd likelihoodCov;

        std::shared_ptr<muq::Modeling::Density> likelihood;

        Eigen::VectorXd priorMean;
        Eigen::MatrixXd priorPrecision;
        std::shared_ptr<muq::Modeling::Density> prior;

        std::shared_ptr<muq::Modeling::DensityProduct> prodDens;

        std::shared_ptr<muq::Modeling::WorkGraph> graph;
        std::shared_ptr<muq::Modeling::ModGraphPiece> joinDens;

        std::shared_ptr<muq::SamplingAlgorithms::SamplingProblem> problem;

        std::shared_ptr<muq::SamplingAlgorithms::SingleChainMCMC> mcmc;

        Eigen::VectorXd get_gWeightsByLinSys()
        {
            List<Eigen::MatrixXd> linSys;
            linSys.resize(2);
            linSys[0] = Theta;
            linSys[1] = deterministicProblem.Tmeas + deterministicProblem.addSol;
            return linSys[0].fullPivLu().solve(linSys[1]);
        }

        Eigen::MatrixXd get_MRFcorrelationMatrix(unsigned TCrows, unsigned TCcols)
        {
            //I assume that the TC vector is organized first by rows and then by cols
            //Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Identity(deterministicProblem.thermocouplesNum, deterministicProblem.thermocouplesNum);
            Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Zero(
                                             deterministicProblem.thermocouplesNum, deterministicProblem.thermocouplesNum);

            for (unsigned row = 0; row < TCrows; row++)
            {
                for (unsigned col = 0; col < TCcols; col++)
                {
                    unsigned vectorID = col + TCrows * row;
                    unsigned Nneigh = 0;

                    if (col != 0)
                    {
                        corrMatrix(vectorID, vectorID - 1) = -1;
                        Nneigh++;
                    }

                    if (col != TCcols - 1)
                    {
                        corrMatrix(vectorID, vectorID + 1) = -1;
                        Nneigh++;
                    }

                    if (row != 0)
                    {
                        corrMatrix(vectorID - TCcols, vectorID) = -1;
                        Nneigh++;
                    }

                    if (row != TCrows - 1)
                    {
                        corrMatrix(vectorID + TCcols, vectorID) = -1;
                        Nneigh++;
                    }

                    corrMatrix(vectorID, vectorID) = Nneigh;
                }
            }

            return corrMatrix;
        }

        void assebleGraph()
        {
            graph->AddNode(weightsPiece, "weights");
            graph->AddNode(ThetaOp, "Theta Product");
            graph->AddNode(likelihood, "likelihood");
            graph->AddNode(prior, "prior");
            graph->AddNode(prodDens, "prodDens");
            graph->AddEdge("weights", 0, "Theta Product", 0);
            graph->AddEdge("Theta Product", 0, "likelihood", 0);
            graph->AddEdge("weights", 0, "prior", 0);
            graph->AddEdge("likelihood", 0, "prodDens", 0);
            graph->AddEdge("prior", 0, "prodDens", 1);
            graph->Visualize("Graph.png");
        }

        std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> Run(
            std::vector<Eigen::VectorXd> start)
        {
            return mcmc->Run(start);
        }
};

#endif
