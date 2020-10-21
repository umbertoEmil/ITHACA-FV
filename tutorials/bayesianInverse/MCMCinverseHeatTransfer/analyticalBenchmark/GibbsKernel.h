#ifndef GIBBSKERNEL_H_
#define GIBBSKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"

#include <iomanip>


namespace muq {
  namespace SamplingAlgorithms {

    /**
      @ingroup MCMCKernels
      @class GibbsKernel
      @brief An implementation of the standard Metropolis-Hastings transition kernel.
      @details <B>Configuration Parameters:</B>

      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Proposal"  | String | - | A string pointing to a block of proposal options. |

     */
    class GibbsKernel : public TransitionKernel {
    public:

      GibbsKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem):
	       TransitionKernel(pt, problem)
          {
          
            // Extract the proposal parts from the ptree
            std::string proposalName = pt.get<std::string>("Proposal");
          
            boost::property_tree::ptree subTree = pt.get_child(proposalName);
            subTree.put("BlockIndex", blockInd);
          
            // Construct the proposal
            proposal = MCMCProposal::Construct(subTree, problem);
            assert(proposal);
          };

      GibbsKernel(boost::property_tree::ptree const& pt,
               std::shared_ptr<AbstractSamplingProblem> problem,
               std::shared_ptr<MCMCProposal> proposalIn) : 
	       TransitionKernel(pt, problem),
	       proposal(proposalIn) {};

      virtual ~GibbsKernel() = default;

      virtual inline std::shared_ptr<MCMCProposal> Proposal() {return proposal;};

      virtual void PostStep(unsigned int const t, std::vector<std::shared_ptr<SamplingState>> const& state)
      {
          proposal->Adapt(t,state);
      };

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t, std::shared_ptr<SamplingState> prevState)
      {
        assert(proposal);
      
        // propose a new point
        std::shared_ptr<SamplingState> prop = proposal->Sample(prevState);
      
        // Aceptance probability
        // accept
        numCalls++;
        if (problem->numBlocksQOI > 0) {
          prop->meta["QOI"] = problem->QOI();
        }
        numAccepts++;
      
        prop->meta["IsProposal"] = false;
        return std::vector<std::shared_ptr<SamplingState>>(1, prop);
      };

      virtual void PrintStatus(std::string prefix) const 
      {
        std::stringstream msg;
        msg << std::setprecision(2);
	if(numCalls > 0)
	{
            msg << prefix << "Acceptance Rate = "  << 100.0*double(numAccepts)/double(numCalls) << "%";
	}
      
        std::cout << msg.str() << std::endl;
      };


      virtual inline double AcceptanceRate() const {return double(numAccepts)/double(numCalls);};


    protected:
      std::shared_ptr<MCMCProposal> proposal;

      unsigned int numCalls = 0;
      unsigned int numAccepts = 0;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
