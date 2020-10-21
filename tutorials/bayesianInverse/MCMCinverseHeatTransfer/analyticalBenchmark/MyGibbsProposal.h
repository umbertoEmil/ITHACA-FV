#ifndef MYGIBBSPROPOSAL_H_
#define MYGIBBSPROPOSAL_H_

#include <iostream>
#include "interpolation.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
#include "inverseLaplacianProblem.H"
#include "reducedInverseLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"

#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"

namespace muq {
  namespace SamplingAlgorithms {

    /** @ingroup MCMCProposals
        @class MyGibbsProposal
        @brief Implementation of the classic Random Walk Metropolis proposal
        @details <B>Configuration Parameters:</B>

        Parameter Key | Type | Default Value | Description |
        ------------- | ------------- | ------------- | ------------- |
        "ProposalVariance"  | Double | - | The variance of an isotropic random walk proposal. |
    */
    class MyGibbsProposal : public MCMCProposal {
    private:
      

    public:
      MyGibbsProposal(boost::property_tree::ptree const& pt,
                 std::shared_ptr<AbstractSamplingProblem> prob):
                       MCMCProposal(pt,prob),
		       numBlocks(pt.get("numBlocks", 1)),
		       sigma(pt.get("sigma", 0.0)),
		       alpha(pt.get("alpha", 0.0))
      {
        unsigned int problemDim = prob->blockSizes(blockInd);
        Theta  = ITHACAstream::readMatrix(folderOffline + "Theta_mat.txt");
        addSol = ITHACAstream::readMatrix(folderOffline + "addSol_mat.txt");
        Tmeas  = ITHACAstream::readMatrix(outputFolder + "Tmeas_mat.txt");
	MRFcorrelationMatrix = Eigen::MatrixXd::Identity(16,16);
      };


      virtual ~MyGibbsProposal() = default;

    protected:

      /// The proposal distribution
      std::shared_ptr<muq::Modeling::GaussianBase> proposal;

      const unsigned numBlocks;
      //Eigen::Ref<const Eigen::VectorXd> Tmeas; 
      const double sigma;
      const double alpha;
      word folderOffline = "./ITHACAoutput/offlineParamBC/";//TODO pass it in the pt tree
      word outputFolder = "./ITHACAoutput/";
      Eigen::MatrixXd Theta;
      Eigen::VectorXd addSol;
      Eigen::VectorXd Tmeas;
      Eigen::MatrixXd MRFcorrelationMatrix;
      
      virtual std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override {
        assert(currentState->state.size()>blockInd);
      
        // the mean of the proposal is the current point
        std::vector<Eigen::VectorXd> props = currentState->state;
	Eigen::VectorXd state(Tmeas.size());
	unsigned it = 0;
	for(unsigned it = 0; it < props.size(); it++)
	{
	    assert(props.at(it).size() == 1);
	    state(it) = props.at(it)(0);
	}


	double a_i = (Theta.col(blockInd).squaredNorm() + 2 * alpha * MRFcorrelationMatrix(blockInd, blockInd)) / (sigma * sigma);
	Eigen::MatrixXd ThetaNoCol = Theta;
	ThetaNoCol.col(blockInd).setZero();
	Eigen::MatrixXd WnoCol = MRFcorrelationMatrix;
	WnoCol.col(blockInd).setZero();
	Eigen::MatrixXd WnoRow = MRFcorrelationMatrix;
	WnoRow.row(blockInd).setZero();

	Eigen::VectorXd mu_s = Tmeas + addSol - ThetaNoCol * state; 
        double mu_p = (WnoRow.col(blockInd).adjoint() * state)(0);
	mu_p +=(WnoCol.row(blockInd) * state)(0);


	double b_i = 2 / (sigma * sigma) * ((Theta.col(blockInd).adjoint() * mu_s)(0) - alpha * mu_p);  

	Eigen::VectorXd mu_i = b_i / (2 * a_i) * Eigen::VectorXd::Ones(1);
	Eigen::VectorXd cov_i = 1 / a_i * Eigen::VectorXd::Ones(1); 
	//std::cout << "Theta.col(blockInd) = " << Theta.col(blockInd) << std::endl;
	//std::cout << "Theta.col(blockInd).squaredNorm() = " << Theta.col(blockInd).squaredNorm() << std::endl << std::endl;
	//std::cout << "a_" << blockInd<<" = " << a_i << std::endl;
	//std::cout << "mu_s" << blockInd<<" = " << mu_s << std::endl;
	//std::cout << "WnoRow.col(blockInd).adjoint() * state = " << WnoRow.col(blockInd).adjoint() * state << std::endl;
	//std::cout << "WnoCol.row(blockInd) * state = " << WnoCol.row(blockInd) * state << std::endl;
	//std::cout << "mu_p" << blockInd<<" = " << mu_p << std::endl;
	//std::cout << "Theta.row(blockInd) * mu_s = " << Theta.row(blockInd) * mu_s << std::endl;
	//std::cout << "alpha * mu_p = " << alpha * mu_p << std::endl;
	//std::cout << "((Theta.row(blockInd) * mu_s)(0) - alpha * mu_p) = " <<
	//    ((Theta.row(blockInd) * mu_s)(0) - alpha * mu_p) << std::endl;
	//std::cout << "sigma = " << sigma << std::endl; 
	//std::cout << "alpha = " << alpha << std::endl; 
	//std::cout << "a_" << blockInd<<" = " << a_i << std::endl;
	//std::cout << "mu_p" << blockInd<<" = " << mu_p << std::endl;
	//std::cout << "b_" << blockInd<<" = " << b_i << std::endl;
	//std::cout << "mu_" << blockInd<<" = " << mu_i << std::endl;
	//std::cout << "cov_" << blockInd<<" = " << cov_i << std::endl << std::endl;

        assert(props.size()>blockInd);
      
        // created a Gaussian with scaled identity covariance
	proposal = std::make_shared<muq::Modeling::Gaussian>(mu_i, cov_i);
        Eigen::VectorXd prop = proposal->Sample();

        props.at(blockInd) = prop;
	
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
