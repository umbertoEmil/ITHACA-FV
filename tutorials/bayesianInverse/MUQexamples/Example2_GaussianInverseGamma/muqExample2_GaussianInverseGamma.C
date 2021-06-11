/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    Example of a heat transfer Reduction Problem
SourceFiles
    14optimization.C
\*---------------------------------------------------------------------------*/

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"


#include  <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;


int main()
{
    auto varPiece = std::make_shared<muq::Modeling::IdentityOperator>(1);
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    auto gaussDens = std::make_shared<muq::Modeling::Gaussian>(mu,
                     muq::Modeling::Gaussian::DiagCovariance)->AsDensity();
    std::cout << "Gaussian piece has " << gaussDens->inputSizes.size()
              << " inputs with sizes " << gaussDens->inputSizes.transpose() << std::endl;
    const double alpha = 2.5;
    const double beta = 1.0;
    auto varDens = std::make_shared<muq::Modeling::InverseGamma>(alpha,
                   beta)->AsDensity();
    auto prodDens = std::make_shared<muq::Modeling::DensityProduct>(2);
    auto replOp = std::make_shared<muq::Modeling::ReplicateOperator>(1, 2);
    auto graph = std::make_shared<muq::Modeling::WorkGraph>();
    graph->AddNode(gaussDens, "Gaussian Density");
    graph->AddNode(varPiece, "Variance");
    graph->AddNode(varDens, "Variance Density");
    graph->AddNode(prodDens, "Joint Density");
    graph->AddNode(replOp, "Replicated Variance");
    graph->AddEdge("Variance", 0, "Replicated Variance", 0);
    graph->AddEdge("Replicated Variance", 0, "Gaussian Density", 1);
    graph->AddEdge("Variance", 0, "Variance Density", 0);
    graph->AddEdge("Gaussian Density", 0, "Joint Density", 0);
    graph->AddEdge("Variance Density", 0, "Joint Density", 1);
    graph->Visualize("DensityGraph.png");
    auto jointDens = graph->CreateModPiece("Joint Density");
    auto problem = std::make_shared<muq::SamplingAlgorithms::SamplingProblem>
                   (jointDens);
    pt::ptree pt;
    pt.put("NumSamples", 1e5); // number of MCMC steps
    pt.put("BurnIn", 1e4);
    pt.put("KernelList",
           "Kernel1,Kernel2"); // Name of block that defines the transition kernel
    pt.put("Kernel1.Method", "MHKernel"); // Name of the transition kernel class
    pt.put("Kernel1.Proposal",
           "MyProposal"); // Name of block defining the proposal distribution
    pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
    pt.put("Kernel1.MyProposal.ProposalVariance",
           0.5); // Variance of the isotropic MH proposal
    pt.put("Kernel2.Method", "MHKernel"); // Name of the transition kernel class
    pt.put("Kernel2.Proposal",
           "GammaProposal"); // Name of block defining the proposal distribution
    pt.put("Kernel2.GammaProposal.Method", "InverseGammaProposal");
    pt.put("Kernel2.GammaProposal.InverseGammaNode", "Variance Density");
    pt.put("Kernel2.GammaProposal.GaussianNode", "Gaussian Density");
    auto mcmc = muq::SamplingAlgorithms::MCMCFactory::CreateSingleChain(pt,
                problem);
    std::vector<Eigen::VectorXd> startPt(2);
    startPt.at(0) = mu; // Start the Gaussian block at the mean
    startPt.at(1) = Eigen::VectorXd::Ones(
                        1); // Set the starting value of the variance to 1
    std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps = mcmc->Run(
                startPt);
    Eigen::VectorXd sampMean = samps->Mean();
    std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
    Eigen::VectorXd sampVar = samps->Variance();
    std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
    Eigen::MatrixXd sampCov = samps->Covariance();
    std::cout << "\nSample Covariance = \n" << sampCov << std::endl;
    Eigen::VectorXd sampMom3 = samps->CentralMoment(3);
    std::cout << "\nSample Third Moment = \n" << sampMom3 << std::endl << std::endl;
    return 0;
}
