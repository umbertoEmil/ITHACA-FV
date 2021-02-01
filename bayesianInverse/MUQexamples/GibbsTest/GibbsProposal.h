#ifndef GIBBSPROPOSAL_H_
#define GIBBSPROPOSAL_H_

#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class GibbsProposal
        @brief Implementation of the classic Random Walk Metropolis proposal
        @details <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "ProposalVariance"  | Double | - | The variance of an isotropic random walk proposal. |
    */
    class GibbsProposal : public MCMCProposal {
    private:
      

    public:

      GibbsProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob):
                       MCMCProposal(pt,prob),
		       numBlocks(pt.get("numBlocks", 1)),
		       rho(pt.get("rho", 1.0))
      {
        std::cout << "numBlocks =" << numBlocks<< std::endl;
        unsigned int problemDim = prob->blockSizes(blockInd);
      
        // compute the (diagonal) covariance for the proposal
        const Eigen::VectorXd cov = pt.get("ProposalVariance", 1.0)*
                                    Eigen::VectorXd::Ones(problemDim);
	
      
        // created a Gaussian with scaled identity covariance
        proposal = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(problemDim), cov);
      };

      GibbsProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob,
                 std::shared_ptr<muq::Modeling::GaussianBase> proposalIn):
		 MCMCProposal(pt,prob), proposal(proposalIn),
		 numBlocks(pt.get("numBlocks", 1)),
		 rho(pt.get("rho", 1.0)) {};

      virtual ~GibbsProposal() = default;

    protected:

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::GaussianBase> proposal;

      const unsigned numBlocks;
      const double rho;
      
      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override {
        assert(currentState->state.size()>blockInd);
      
        // the mean of the proposal is the current point
        std::vector<Eigen::VectorXd> props = currentState->state;

        assert(props.size()>blockInd);
	Eigen::VectorXd xc;
	if(blockInd < numBlocks - 1)
	{
	    xc = currentState->state.at(blockInd + 1);
	}
	else
	{
	    xc = currentState->state.at(0);
	}
      
        Eigen::VectorXd prop = proposal->Sample();
        props.at(blockInd) = xc * rho + std::sqrt(1 - rho * rho) * prop;
	
	//std::cout << "debug block " << blockInd << std::endl;
	//for (const auto& i: props)
	//{
	//    std::cout << i << std::endl;
	//}
	//    std::cout << std::endl;
        // store the new state in the output
        return std::make_shared<SamplingState>(props, 1.0);
      };


      virtual double LogDensity(std::shared_ptr<SamplingState> const& currState,
                                std::shared_ptr<SamplingState> const& propState) override {

        Eigen::VectorXd diff = propState->state.at(blockInd)-currState->state.at(blockInd);
        return proposal->LogDensity(diff);//, std::pair<boost::any, Gaussian::Mode>(conditioned->state.at(blockInd), Gaussian::Mode::Mean));
      };


    };

  } // namespace SamplingAlgoirthms
} // namespace muq


#endif
