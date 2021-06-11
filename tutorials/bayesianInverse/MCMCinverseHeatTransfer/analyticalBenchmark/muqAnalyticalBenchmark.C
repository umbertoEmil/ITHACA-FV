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

#include "MyGibbsProposal.h"
#include "GibbsKernel.h"

#include  <boost/property_tree/ptree.hpp>

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
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Foam2Eigen.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "analyticalBenchmark.H"
#include "bayesianAnalyticalBenchmark.H"



namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

class appendVectors : public muq::Modeling::ModPiece
{
    public:
        appendVectors(int inputSize) : muq::Modeling::ModPiece(1 *
                    Eigen::VectorXi::Ones(inputSize),    // inputSizes = [1,1]
                    inputSize * Eigen::VectorXi::Ones(1)) {}; // outputSizes = [2]
    protected:
        virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const&
                                  inputs) override
        {
            Eigen::VectorXd xjoin;

            for (unsigned inputI = 0; inputI < inputs.size(); inputI++)
            {
                Eigen::VectorXd input = inputs.at(inputI);
                xjoin.resize(xjoin.size() +  input.size());
                xjoin << xjoin, input;
            }

            // Resize the outputs vector (which lives in the ModPiece base class)
            outputs.resize(1);
            outputs.at(0) = xjoin;
        }
};

Eigen::VectorXd expList(double first, double last, Eigen::DenseIndex n)
{
    Eigen::VectorXd vector(n);
    double m = (double) 1 / (n - 1);
    double quotient = std::pow(last / first, m);
    vector(0) = first;

    for (Eigen::DenseIndex i = 1; i < n;
            i++) // DenseIndex is just a typedef ptrdiff_t from the Eigen library
    {
        vector(i) = vector(i - 1) * quotient;
    }

    return vector;
}


int main(int argc, char* argv[])
{
    solverPerformance::debug = 1; //No verbose output
    double time;
    analyticalBenchmark example(argc, argv);
    //Setting parameters for the analytical benchmark
    word outputFolder = "./ITHACAoutput/";
    word resultsFolder = "./ITHACAoutput/muqResults/";
    word resultsFolder_determ = "./ITHACAoutput/deterministic/";
    // Reading tests to perform
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    // Reading parameters from ITHACAdict
    example.a = para->ITHACAdict->lookupOrDefault<double>("a", 0.0);
    example.b = para->ITHACAdict->lookupOrDefault<double>("b", 0.0);
    example.c = para->ITHACAdict->lookupOrDefault<double>("c", 0.0);
    example.d = para->ITHACAdict->lookupOrDefault<double>("d", 0.0);
    word solver = para->ITHACAdict->lookupOrDefault<word>("linSysSolver",
                  "fullPivLU");
    label TSVDtruncation =
        para->ITHACAdict->lookupOrDefault<label>("TSVDtruncation", 3);
    label TCrows = para->ITHACAdict->lookupOrDefault<label>("TCrows", 0);
    label TCcols = para->ITHACAdict->lookupOrDefault<label>("TCcols", 0);
    M_Assert(TCrows > 0
             && TCcols > 0,
             "Set up TCrows and TCcols with the rows and cols of the thermocouples");
    const unsigned int Nsamples =
        para->ITHACAdict->lookupOrDefault<int>("MCMCsamples", 0);
    const unsigned int BurnIn = para->ITHACAdict->lookupOrDefault<int>("MCMCBurnIn",
                                0);
    const unsigned int printLevel =
        para->ITHACAdict->lookupOrDefault<int>("printLevel", 0);
    double measurementsStdDevPerCent =
        para->ITHACAdict->lookupOrDefault<double>("measurementsStdDev", 0);
    double regularizationParameter =
        para->ITHACAdict->lookupOrDefault<double>("regularizationParameter", 0);
    double ProposalVariance =
        para->ITHACAdict->lookupOrDefault<double>("ProposalVariance", 0);
    example.k = para->ITHACAdict->lookupOrDefault<double>("thermalConductivity", 0);
    M_Assert(example.k > 0, "thermalConductivity, k, not specified");
    example.H = para->ITHACAdict->lookupOrDefault<double>("heatTranferCoeff", 0);
    M_Assert(example.H > 0, "Heat transfer coeff, H, not specified");
    label Ntests = para->ITHACAdict->lookupOrDefault<double>("NumberErrorTests",
                   100);
    double Tf = para->ITHACAdict->lookupOrDefault<double>("Tf", 300.0);
    double refGrad = para->ITHACAdict->lookupOrDefault<double>("refGrad", 0.0);
    double valueFraction =
        para->ITHACAdict->lookupOrDefault<double>("valueFraction",
                0.0);
    //Tests to perform
    double alphaMin =
        para->ITHACAdict->lookupOrDefault<double>("regularizationParameterMin", 0);
    double alphaMax =
        para->ITHACAdict->lookupOrDefault<double>("regularizationParameterMax", 0);
    int alphaSize =
        para->ITHACAdict->lookupOrDefault<int>("regularizationParameterSamples", 0);
    bool GibbsSingleRun = para->ITHACAdict->lookupOrDefault<bool>("GibbsSingleRun",
                          0);
    bool MHMCMCRun = para->ITHACAdict->lookupOrDefault<bool>("MHMCMCRun", 0);
    bool MHMCMCRunOld = para->ITHACAdict->lookupOrDefault<bool>("MHMCMCRunOld", 0);
    bool proposalVarianceTest =
        para->ITHACAdict->lookupOrDefault<bool>("proposalVarianceTest", 0);
    bool regularizationParameterTest =
        para->ITHACAdict->lookupOrDefault<bool>("regularizationParameterTest", 0);
    // setting analytical solution
    volScalarField T_true(example._T());

    for (label i = 0; i < T_true.internalField().size(); i++)
    {
        auto cx = T_true.mesh().C()[i].component(vector::X);
        auto cy = T_true.mesh().C()[i].component(vector::Y);
        auto cz = T_true.mesh().C()[i].component(vector::Z);
        T_true.ref()[i] = example.a * cx * cx + example.b * cx * cy + example.c * cy -
                          example.a * cz * cz + example.c;
    }

    fvMesh& mesh = example._mesh();
    example.hotSide_ind = mesh.boundaryMesh().findPatchID("hotSide");
    label hotSideSize = T_true.boundaryField()[example.hotSide_ind].size();
    example.g.resize(hotSideSize);
    forAll(example.g, faceI)
    {
        scalar faceX =
            mesh.boundaryMesh()[example.hotSide_ind].faceCentres()[faceI].x();
        scalar faceZ =
            mesh.boundaryMesh()[example.hotSide_ind].faceCentres()[faceI].z();
        example.g[faceI] = example.k * (example.b * faceX);// + example.c * faceZ) ;
    }
    example.gTrue = example.g; //-example.heatFlux_hotSide / example.k;
    volScalarField gTrueField = example.list2Field(example.gTrue);
    ITHACAstream::exportSolution(gTrueField,
                                 "1", resultsFolder,
                                 "gTrue");
    example.solveTrue();
    volScalarField& T(example._T());
    Info << "Exporting analytical solution" << endl;
    ITHACAstream::exportSolution(T_true, "1", "./ITHACAoutput/true/",
                                 "analyticalSol");
    // Setting up the thermocouples
    example.readThermocouples();
    example.Tmeas = example.fieldValueAtThermocouples(T_true);
    ITHACAstream::exportMatrix(example.Tmeas, "Tmeas", "eigen", outputFolder);
    double measurementsStdDev = measurementsStdDevPerCent * example.Tmeas.mean();
    // Introducing error in the measurements
    Eigen::VectorXd TmeasTrue = example.Tmeas;
    auto measErrorDist = std::make_shared<muq::Modeling::Gaussian>
                         (Eigen::VectorXd::Zero(TmeasTrue.size()),
                          measurementsStdDev * measurementsStdDev * Eigen::MatrixXd::Identity(
                              TmeasTrue.size(), TmeasTrue.size()));
    example.Tmeas += measErrorDist->Sample();
    Info << endl;
    Info << "*********************************************************" << endl;
    Info << "Computing offline matrices and vectors" << endl;
    Info << endl;
    example.set_gParametrized("rbf", 0.7);
    example.parameterizedBCoffline();
    M_Assert(TCrows * TCcols == example.Tmeas.size(),
             "Control number of thermocouples rows and columns");
    //Eigen::MatrixXd MRFcorrelationMatrix = example.get_MRFcorrelationMatrix(TCrows, TCcols);
    //ITHACAstream::exportMatrix(MRFcorrelationMatrix, "MRFcorrelationMatrix",
    //                                   "eigen", outputFolder);
    Info << "*********************************************************" << endl;
    Info << endl;
    Info << "*********************************************************" << endl;
    Info << "Deterministic approach" << endl;
    Info << "Solver " << solver << endl;
    Info << endl;
    example.parameterizedBC(solver, TSVDtruncation);
    volScalarField gParametrizedField = example.list2Field(example.g);
    ITHACAstream::exportSolution(gParametrizedField,
                                 std::to_string(1),
                                 resultsFolder_determ,
                                 "gParametrized");
    volScalarField& Tdet(example._T());
    ITHACAstream::exportSolution(Tdet,
                                 std::to_string(1),
                                 resultsFolder_determ,
                                 "T_param");
    Info << "residualNorm = " << Foam::sqrt(
             example.residual.squaredNorm()) << endl;
    Eigen::MatrixXd A = example.Theta.transpose() * example.Theta;
    example.postProcess(resultsFolder_determ, "gParametrized");
    Info << "debug: gWeights = " << example.gWeights << endl;
    Info << "*********************************************************" << endl;
    Info << endl;
    Info << "*********************************************************" << endl;
    Info << "Bayesian Approach" << endl;
    REGISTER_TRANSITION_KERNEL(GibbsKernel)
    REGISTER_MCMC_PROPOSAL(MyGibbsProposal)

    if (MHMCMCRun)
    {
        word MHMCMCresultsFolder = "./ITHACAoutput/MHMCMCresults/";
        //unsigned int AdaptSteps = 1000;
        //unsigned int AdaptStart = 2000;
        //double AdaptScale = 0.5;
        //// parameters for the sampler
        pt::ptree pt;
        pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
        Info << "debug Nsamples = " << Nsamples << endl;
        pt.put("BurnIn", BurnIn);
        pt.put("PrintLevel", printLevel);
        pt.put("KernelList",
               "Kernel1"); // Name of block that defines the transition kernel
        pt.put("Kernel1.Method", "MHKernel"); // Name of the transition kernel class
        pt.put("Kernel1.Proposal",
               "MyProposal"); // Name of block defining the proposal distribution
        pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
        pt.put("Kernel1.MyProposal.ProposalVariance",
               ProposalVariance); // Variance of the isotropic MH proposal
        //pt.put("Kernel1.MyProposal.AdaptSteps", AdaptSteps); // Variance of the isotropic MH proposal
        //pt.put("Kernel1.MyProposal.AdaptStart", AdaptStart);
        //pt.put("Kernel1.MyProposal.AdaptScale", AdaptScale);
        bayesianAnalyticalBenchmark bayesianExample(example, measurementsStdDev,
                regularizationParameter, TCrows, TCcols, pt);
        // starting point
        std::vector<Eigen::VectorXd> start(1);
        start.at(0) = Eigen::VectorXd::Zero(bayesianExample.Nweights);
        std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps =
            bayesianExample.Run(start);
        Info << "Post processing \n" << endl;
        Eigen::VectorXd sampMean = samps->Mean();
        std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
        std::cout << "Theta * Sample Mean = \n" << (example.Theta*
                  sampMean).transpose() << std::endl;
        std::cout << "example.addSol + example.Tmeas = \n" << (example.addSol +
                  example.Tmeas).transpose() << std::endl;
        Eigen::VectorXd sampVar = samps->Variance();
        std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
        //Eigen::VectorXd ESS = samps->ESS();
        //std::cout << "\nEffective sample size = \n" << ESS.transpose() << std::endl;
        //MAP
        Eigen::MatrixXd logTarget = samps->GetMeta("LogTarget");
        samps->WriteToFile("MHresults");
        example.parameterizedBCpostProcess(sampMean);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(1),
                                     MHMCMCresultsFolder,
                                     "g_MH");
        volScalarField& T_MH(example._T());
        ITHACAstream::exportSolution(T_MH,
                                     std::to_string(1),
                                     MHMCMCresultsFolder,
                                     "T_MH");
        example.postProcess(MHMCMCresultsFolder, "g_MH");
    }

    if (proposalVarianceTest)
    {
        Eigen::VectorXd alpha = Eigen::VectorXd::LinSpaced(alphaSize, alphaMin,
                                alphaMax);
        Eigen::VectorXd Js(alphaSize);
        word MHMCMCregParResultsFolder = "./ITHACAoutput/MHMCMCpropVarTestResults/";

        for (int alphaI = 0; alphaI < alphaSize; alphaI++)
        {
            //// parameters for the sampler
            pt::ptree pt;
            pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
            pt.put("BurnIn", BurnIn);
            pt.put("PrintLevel", printLevel);
            pt.put("KernelList",
                   "Kernel1"); // Name of block that defines the transition kernel
            pt.put("Kernel1.Method", "MHKernel"); // Name of the transition kernel class
            pt.put("Kernel1.Proposal",
                   "MyProposal"); // Name of block defining the proposal distribution
            pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
            pt.put("Kernel1.MyProposal.ProposalVariance",
                   alpha(alphaI)); // Variance of the isotropic MH proposal
            bayesianAnalyticalBenchmark bayesianExample(example, measurementsStdDev,
                    regularizationParameter, TCrows, TCcols, pt);
            // starting point
            std::vector<Eigen::VectorXd> start(1);
            start.at(0) = Eigen::VectorXd::Zero(bayesianExample.Nweights);
            std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps =
                bayesianExample.Run(start);
            Info << "Post processing \n" << endl;
            Eigen::VectorXd sampMean = samps->Mean();
            std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
            Eigen::VectorXd sampVar = samps->Variance();
            std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
            Eigen::VectorXd ESS = samps->ESS();
            std::cout << "\nEffective sample size = \n" << ESS.transpose() << std::endl;
            example.parameterizedBCpostProcess(sampMean);
            gParametrizedField = example.list2Field(example.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(alphaI + 1),
                                         MHMCMCregParResultsFolder,
                                         "g_MH");
            volScalarField& T_MH(example._T());
            ITHACAstream::exportSolution(T_MH,
                                         std::to_string(alphaI + 1),
                                         MHMCMCregParResultsFolder,
                                         "T_MH");
            Js(alphaI) = example.J;
        }

        example.postProcess(MHMCMCregParResultsFolder, "g_MH");
        Eigen::MatrixXd Jmatrix(Js.size(), 2);
        Jmatrix << alpha, Js;
        ITHACAstream::exportMatrix(Jmatrix, "Js", "eigen", MHMCMCregParResultsFolder);
    }

    if (regularizationParameterTest)
    {
        Eigen::VectorXd alpha = expList(alphaMin, alphaMax, alphaSize);
        Eigen::VectorXd Js(alphaSize);
        word MHMCMCregParResultsFolder = "./ITHACAoutput/MHMCMCregParTestResults/";
        //// parameters for the sampler
        pt::ptree pt;
        pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
        pt.put("BurnIn", BurnIn);
        pt.put("PrintLevel", printLevel);
        pt.put("KernelList",
               "Kernel1"); // Name of block that defines the transition kernel
        pt.put("Kernel1.Method", "MHKernel"); // Name of the transition kernel class
        pt.put("Kernel1.Proposal",
               "MyProposal"); // Name of block defining the proposal distribution
        pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
        pt.put("Kernel1.MyProposal.ProposalVariance",
               ProposalVariance); // Variance of the isotropic MH proposal

        for (int alphaI = 0; alphaI < alphaSize; alphaI++)
        {
            Info << "\n\nIteration " << alphaI + 1 << ", alpha = " << alpha(alphaI) << endl;
            auto bayesianExample = std::make_shared<bayesianAnalyticalBenchmark>(example,
                                   measurementsStdDev, alpha(alphaI), TCrows, TCcols, pt);
            // starting point
            std::vector<Eigen::VectorXd> start(1);
            start.at(0) = Eigen::VectorXd::Zero(bayesianExample->Nweights);
            std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps =
                bayesianExample->Run(start);
            Info << "Post processing \n" << endl;
            Eigen::VectorXd sampMean = samps->Mean();
            std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
            Eigen::VectorXd sampVar = samps->Variance();
            std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
            example.parameterizedBCpostProcess(sampMean);
            gParametrizedField = example.list2Field(example.g);
            ITHACAstream::exportSolution(gParametrizedField,
                                         std::to_string(alphaI + 1),
                                         MHMCMCregParResultsFolder,
                                         "g_MH");
            volScalarField& T_MH(example._T());
            ITHACAstream::exportSolution(T_MH,
                                         std::to_string(alphaI + 1),
                                         MHMCMCregParResultsFolder,
                                         "T_MH");
            Js(alphaI) = example.J;
            samps.reset();
        }

        //example.postProcess(MHMCMCregParResultsFolder, "g_MH");
        Eigen::MatrixXd Jmatrix(Js.size(), 2);
        Jmatrix << alpha, Js;
        ITHACAstream::exportMatrix(Jmatrix, "Js", "eigen", MHMCMCregParResultsFolder);
    }

    //    if(regularizationParameterTest)
    //    {
    //        Eigen::VectorXd alpha = Eigen::VectorXd::LinSpaced(alphaSize, alphaMin, alphaMax);
    //  Eigen::VectorXd Js(alphaSize);
    //        for(int alphaI = 0; alphaI < alphaSize; alphaI++)
    //        {
    //
    //            // parameters for the sampler
    //            pt::ptree pt;
    //            pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
    //            Info << "debug Nsamples = " << Nsamples << endl;
    //            pt.put("BurnIn", BurnIn);
    //            pt.put("PrintLevel",printLevel);
    //            std::string kernelsList;
    //            unsigned int numBlocks = example.gWeights.size();
    //            forAll(example.gWeights, kernelI)
    //            {
    //                kernelsList.append("Kernel" + std::to_string(kernelI));
    //                if(kernelI < example.gWeights.size() - 1)
    //                {
    //                    kernelsList.append(",");
    //                }
    //                pt.put("Kernel" + std::to_string(kernelI) + ".Method","GibbsKernel");
    //                pt.put("Kernel" + std::to_string(kernelI) + ".Proposal", "MyProposal");
    //                pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.Method", "MyGibbsProposal");
    //                pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.numBlocks", numBlocks);
    //                pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.sigma", measurementsStdDev);
    //                pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.alpha", alpha(alphaI));
    //            }
    //            pt.put("KernelList", kernelsList); // the transition kernel
    //
    //
    //            auto graph = std::make_shared<muq::Modeling::WorkGraph>();
    //            auto appendVect = std::make_shared<appendVectors>(numBlocks);
    //
    //            graph->AddNode(appendVect, "appendVectors");
    //            graph->Visualize("Graph.png");
    //
    //            auto dens = graph->CreateModPiece("appendVectors");
    //            // create a sampling problem
    //            auto problem = std::make_shared<muq::SamplingAlgorithms::SamplingProblem>(dens);
    //
    //            // starting point
    //            std::vector<Eigen::VectorXd> start(numBlocks);
    //            forAll(example.gWeights, weightI)
    //            {
    //                 start.at(weightI) = Eigen::VectorXd::Zero(1);
    //            }
    //
    //            // evaluate
    //            // create an instance of MCMC
    //            auto mcmc = std::make_shared<muq::SamplingAlgorithms::MySingleChainMCMC>(pt,problem);
    //
    //            std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps = mcmc->Run(start);
    //
    //            Info << "Post processing \n" << endl;
    //            Eigen::VectorXd sampMean = samps->Mean();
    //            std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
    //
    //            Eigen::VectorXd sampVar = samps->Variance();
    //            std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
    //
    //      //samps->WriteToFile("results" + std::to_string(alphaI));
    //
    //            example.parameterizedBCpostProcess(sampMean);
    //            gParametrizedField = example.list2Field(example.g);
    //            ITHACAstream::exportSolution(gParametrizedField,
    //                                         std::to_string(alphaI + 1),
    //                                         resultsFolder,
    //                                         "g_Gibbs");
    //            volScalarField& T_Gibbs(example._T());
    //            ITHACAstream::exportSolution(T_Gibbs,
    //                                         std::to_string(alphaI + 1),
    //                                         resultsFolder,
    //                                         "T_Gibbs");
    //      Js(alphaI) = example.J;
    //        }
    //        example.postProcess(resultsFolder, "g_Gibbs");
    //  Eigen::MatrixXd Jmatrix(Js.size(),2);
    //  Jmatrix << alpha, Js;
    //        ITHACAstream::exportMatrix(Jmatrix, "Js", "eigen", resultsFolder);
    //    }
    //
    //
    if (GibbsSingleRun)
    {
        word GibbsResultsFolder = "./ITHACAoutput/GibbsResults/";
        // parameters for the sampler
        pt::ptree pt;
        pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
        Info << "debug Nsamples = " << Nsamples << endl;
        pt.put("BurnIn", BurnIn);
        pt.put("PrintLevel", printLevel);
        std::string kernelsList;
        unsigned int numBlocks = example.gWeights.size();
        forAll(example.gWeights, kernelI)
        {
            kernelsList.append("Kernel" + std::to_string(kernelI));

            if (kernelI < example.gWeights.size() - 1)
            {
                kernelsList.append(",");
            }

            pt.put("Kernel" + std::to_string(kernelI) + ".Method", "GibbsKernel");
            pt.put("Kernel" + std::to_string(kernelI) + ".Proposal", "MyProposal");
            pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.Method",
                   "MyGibbsProposal");
            pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.numBlocks", numBlocks);
            pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.sigma",
                   measurementsStdDev);
            pt.put("Kernel" + std::to_string(kernelI) + ".MyProposal.alpha",
                   regularizationParameter);
        }
        pt.put("KernelList", kernelsList); // the transition kernel
        auto graph = std::make_shared<muq::Modeling::WorkGraph>();
        auto appendVect = std::make_shared<appendVectors>(numBlocks);
        graph->AddNode(appendVect, "appendVectors");
        graph->Visualize("Graph.png");
        auto dens = graph->CreateModPiece("appendVectors");
        // create a sampling problem
        auto problem = std::make_shared<muq::SamplingAlgorithms::SamplingProblem>(dens);
        // starting point
        std::vector<Eigen::VectorXd> start(numBlocks);
        forAll(example.gWeights, weightI)
        {
            start.at(weightI) = Eigen::VectorXd::Zero(1);
        }
        // evaluate
        // create an instance of MCMC
        auto mcmc = std::make_shared<muq::SamplingAlgorithms::MySingleChainMCMC>(pt,
                    problem);
        std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps = mcmc->Run(
                    start);
        Info << "Post processing \n" << endl;
        Eigen::VectorXd sampMean = samps->Mean();
        std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
        Eigen::VectorXd sampVar = samps->Variance();
        std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
        //Eigen::VectorXd ESS = samps->ESS();
        //std::cout << "\nEffective sample size = \n" << ESS.transpose() << std::endl;
        //Eigen::MatrixXd sampCov = samps->Covariance();
        //std::cout << "\nSample Covariance = \n" << sampCov << std::endl;
        //Eigen::VectorXd sampMom3 = samps->CentralMoment(3);
        //std::cout << "\nSample Third Moment = \n" << sampMom3 << std::endl << std::endl;
        samps->WriteToFile("GibbsResults");
        example.parameterizedBCpostProcess(sampMean);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(1),
                                     GibbsResultsFolder,
                                     "g_Gibbs");
        volScalarField& T_Gibbs(example._T());
        ITHACAstream::exportSolution(T_Gibbs,
                                     std::to_string(1),
                                     resultsFolder,
                                     "T_Gibbs");
        example.postProcess(GibbsResultsFolder, "g_Gibbs");
        Info << "debug: gWeights = " << example.gWeights << endl;
    }

    if (MHMCMCRunOld)
    {
        word MHMCMCresultsFolder = "./ITHACAoutput/MHMCMCresultsOld/";
        Eigen::VectorXd weights(example.gWeights.size());
        //// parameters for the sampler
        pt::ptree pt;
        pt.put("NumSamples", Nsamples); // number of Monte Carlo samples
        Info << "debug Nsamples = " << Nsamples << endl;
        pt.put("BurnIn", BurnIn);
        pt.put("PrintLevel", printLevel);
        pt.put("KernelList",
               "Kernel1"); // Name of block that defines the transition kernel
        pt.put("Kernel1.Method", "MHKernel"); // Name of the transition kernel class
        pt.put("Kernel1.Proposal",
               "MyProposal"); // Name of block defining the proposal distribution
        pt.put("Kernel1.MyProposal.Method", "MHProposal"); // Name of proposal class
        pt.put("Kernel1.MyProposal.ProposalVariance",
               ProposalVariance); // Variance of the isotropic MH proposal
        bayesianAnalyticalBenchmark bayesianExample(example, measurementsStdDev,
                regularizationParameter, TCrows, TCcols, pt);
        auto weightsPiece = std::make_shared<IdentityOperator>(weights.size());
        Eigen::MatrixXd Theta = example.Theta;
        std::shared_ptr<muq::Modeling::LinearOperator> ThetaOp =
            muq::Modeling::LinearOperator::Create(Theta);
        //std::cout << "thetaOp piece has " << ThetaOp->inputSizes.size()
        //          << " inputs with sizes " << ThetaOp->inputSizes.transpose() << std::endl;
        double likelihoodStdDev = measurementsStdDev;
        Eigen::MatrixXd likelihoodCov = likelihoodStdDev * likelihoodStdDev *
                                        Eigen::MatrixXd::Identity(example.Tmeas.size(), example.Tmeas.size());
        auto likelihood = std::make_shared<muq::Modeling::Gaussian>
                          (example.addSol + example.Tmeas, likelihoodCov)->AsDensity();
        Eigen::VectorXd priorMean = Eigen::VectorXd::Zero(example.gWeights.size());
        Eigen::MatrixXd priorPrecision = bayesianExample.MRFcorrelationMatrix *
                                         (2.0 * regularizationParameter / (likelihoodStdDev * likelihoodStdDev));
        auto prior = std::make_shared<muq::Modeling::Gaussian>(priorMean,
                     priorPrecision, muq::Modeling::Gaussian::Mode::Precision )->AsDensity();
        auto prodDens = std::make_shared<muq::Modeling::DensityProduct>(2);
        auto graph = std::make_shared<muq::Modeling::WorkGraph>();
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
        auto jointDens = graph->CreateModPiece("prodDens");
        auto problem = std::make_shared<SamplingProblem>(jointDens);
        auto mcmc = std::make_shared<SingleChainMCMC>(pt, problem);
        // starting point
        std::vector<Eigen::VectorXd> start(1);
        start.at(0) = Eigen::VectorXd::Zero(weights.size());
        std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samps = mcmc->Run(
                    start);
        Info << "Post processing \n" << endl;
        Eigen::VectorXd sampMean = samps->Mean();
        std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
        Eigen::VectorXd sampVar = samps->Variance();
        std::cout << "\nSample Variance = \n" << sampVar.transpose() << std::endl;
        //Eigen::VectorXd ESS = samps->ESS();
        //std::cout << "\nEffective sample size = \n" << ESS.transpose() << std::endl;
        samps->WriteToFile("MHresults");
        example.parameterizedBCpostProcess(sampMean);
        gParametrizedField = example.list2Field(example.g);
        ITHACAstream::exportSolution(gParametrizedField,
                                     std::to_string(1),
                                     MHMCMCresultsFolder,
                                     "g_MH");
        volScalarField& T_Gibbs(example._T());
        ITHACAstream::exportSolution(T_Gibbs,
                                     std::to_string(1),
                                     MHMCMCresultsFolder,
                                     "T_MH");
        example.postProcess(MHMCMCresultsFolder, "g_MH");
        Info << "debug: gWeights = " << example.gWeights << endl;
    }

    return 0;
}
