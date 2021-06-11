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
    Example of implementation of Gibbs sampling
SourceFiles
    muqGibbsTest.C 
\*---------------------------------------------------------------------------*/

#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/MALAProposal.h"
#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/DRKernel.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"

#include "GibbsProposal.h"
#include "GibbsKernel.h"

#include  <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

class appendVectors : public muq::Modeling::ModPiece
{
public:
  appendVectors() : muq::Modeling::ModPiece(1*Eigen::VectorXi::Ones(2),    // inputSizes = [1,1]
                                                   2*Eigen::VectorXi::Ones(1)){}; // outputSizes = [2]
protected:
  virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override
  {
    Eigen::VectorXd const& x1 = inputs.at(0);
    Eigen::VectorXd const& x2 = inputs.at(1);
    Eigen::VectorXd xjoin(x1.size() +  x2.size());
    xjoin << x1, x2;
    // Resize the outputs vector (which lives in the ModPiece base class)
    outputs.resize(1);
    outputs.at(0) = xjoin;
  }
};


int main(){
  const unsigned int N = 1e4;

  REGISTER_TRANSITION_KERNEL(GibbsKernel)
  REGISTER_MCMC_PROPOSAL(GibbsProposal)
  // parameters for the sampler
  pt::ptree pt;
  pt.put("NumSamples", N); // number of Monte Carlo samples
  pt.put("PrintLevel",3);
  pt.put("KernelList", "Kernel1,Kernel2"); // the transition kernel
  //pt.put("KernelList", "Kernel1"); // the transition kernel

  std::string kernelString = pt.get<std::string>("KernelList");
  std::vector<std::string> kernelNames = muq::Utilities::StringUtilities::Split(kernelString, ',');

  std::vector<std::shared_ptr<TransitionKernel>> kernelVec;
  unsigned int numBlocks = kernelNames.size();
  const double rho = 0.9;

  // MH Kernel for first block
  pt.put("Kernel1.Method","GibbsKernel");
  pt.put("Kernel1.Proposal", "MyProposal");
  pt.put("Kernel1.MyProposal.Method", "GibbsProposal");
  pt.put("Kernel1.MyProposal.ProposalVariance", 1);
  pt.put("Kernel1.MyProposal.numBlocks", numBlocks);
  pt.put("Kernel1.MyProposal.rho", rho);
  
  pt.put("Kernel2.Method","GibbsKernel");
  pt.put("Kernel2.Proposal", "MyProposal");
  pt.put("Kernel2.MyProposal.Method", "GibbsProposal");
  pt.put("Kernel2.MyProposal.ProposalVariance", 1);
  pt.put("Kernel2.MyProposal.numBlocks", numBlocks);
  pt.put("Kernel2.MyProposal.rho", rho);

  const Eigen::VectorXd mu1 = 0*Eigen::VectorXd::Ones(1);
  const Eigen::VectorXd mu2 = 1*Eigen::VectorXd::Ones(1);

  auto graph = std::make_shared<muq::Modeling::WorkGraph>();
  auto appendVect = std::make_shared<appendVectors>();
  
  graph->AddNode(appendVect, "appendVectors");
  graph->Visualize("UnconnectedGraph.png");

  graph->Visualize("DensityGraph.png");

  auto dens = graph->CreateModPiece("appendVectors");
  // create a sampling problem
  auto problem = std::make_shared<muq::SamplingAlgorithms::SamplingProblem>(dens);

  // starting point
  std::vector<Eigen::VectorXd> start(2);
  start.at(0) = mu1;
  start.at(1) = mu2;

  // evaluate
  // create an instance of MCMC
  auto mcmc = std::make_shared<muq::SamplingAlgorithms::SingleChainMCMC>(pt,problem);

  std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps = mcmc->Run(start);

  
  Eigen::VectorXd sampMean = samps->Mean();
  std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;

  Eigen::VectorXd sampVar = samps->Variance();
  std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;

  Eigen::MatrixXd sampCov = samps->Covariance();
  std::cout << "\nSample Covariance = \n" << sampCov << std::endl;

  Eigen::VectorXd sampMom3 = samps->CentralMoment(3);
  std::cout << "\nSample Third Moment = \n" << sampMom3 << std::endl << std::endl;
  samps->WriteToFile("results");

  return 0;
}
