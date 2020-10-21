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
    Test of the Ensemble Kalman filter implementation
SourceFiles
    EnKFtest2.C 
\*---------------------------------------------------------------------------*/


#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/Utilities/AnyHelpers.h"
#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"

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
#include "muq2ithaca.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"


int main(int argc, char* argv[])
{
    word outputFolder = "./ITHACAoutput/";

    int Nseeds = 10000;

    Eigen::MatrixXd A = ITHACAstream::readMatrix("A_mat.txt");
    Eigen::MatrixXd H = ITHACAstream::readMatrix("observation_mat.txt");
    M_Assert(ITHACAstream::readMatrix("initialState_mat.txt").cols() == 1, "Wrong initialState input");
    Eigen::VectorXd x0 = ITHACAstream::readMatrix("initialState_mat.txt").col(0);
    int stateSize = A.rows();
    int obsSize = H.rows();
    Eigen::VectorXd b  (stateSize);
    b << 1, 2, 3, 4;
    std::cout << "In this tutorial we have a dynamical system in the form:\nx_{n+1} = A * x_{n} + b" << std::endl;
    std::cout << "with A = \n" << A << std::endl;
    std::cout << "and b = \n(" << b.transpose() << ") * sin(t)" << std::endl;

    std::cout << "We observe the state x by mean of the observation matrix \nH = \n" << H << std::endl; 
    std::cout << "The objective is to reconstruct the vector b knowing H, x_0, and A" << std::endl; 
 

    int Ntimes = 51;
    int sampleDeltaStep = 1;
    double endTime = 10;
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(Ntimes, 0, endTime);

    Eigen::VectorXd xOld = x0;
    Eigen::MatrixXd X(stateSize,Ntimes);
    Eigen::MatrixXd forcer(stateSize,Ntimes);
    X.col(0) = x0;
    forcer.col(0) = b *0.0;

    int sampleFlag = sampleDeltaStep;
    int Nsamples = (Ntimes - 1) / sampleDeltaStep;
    int sampleI = 0;
    Eigen::MatrixXd obs(obsSize, Nsamples);

    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        forcer.col(timeI + 1) = b * std::sin(time(timeI+1));
	Eigen::VectorXd xNew = A * xOld + forcer.col(timeI + 1);
	xOld = xNew;
	Eigen::VectorXd dNew = H * xNew;
        X.col(timeI + 1) = xNew;
        sampleFlag--;
        if(sampleFlag == 0)
        {
            sampleFlag = sampleDeltaStep;
            obs.col(sampleI) = dNew;
            sampleI++;
        }
    }
    M_Assert(Nsamples == sampleI, "Something went wrong in the sampling");

    ITHACAstream::exportMatrix(time, "time", "eigen", outputFolder);
    ITHACAstream::exportMatrix(X, "X", "eigen", outputFolder);
    ITHACAstream::exportMatrix(forcer, "forcer", "eigen", outputFolder);

    Eigen::VectorXd x = x0;


    Eigen::MatrixXd meas_cov = Eigen::MatrixXd::Identity(obsSize, obsSize) * 0.5;
    auto measNoise = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(obsSize), meas_cov); 

    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    posteriorMean.col(0) *= 0.0;
    Eigen::MatrixXd stateRec = posteriorMean;
    Eigen::MatrixXd minConfidence = posteriorMean;
    minConfidence.col(0) *= 0.0;
    Eigen::MatrixXd maxConfidence = minConfidence;
    stateRec.col(0) = x0;
    //Eigen::VectorXd prior_mu = Eigen::VectorXd::Zero(stateSize);//posteriorMean.col(timeI);

    Eigen::MatrixXd singleVariableSamples(Ntimes, Nseeds);
    singleVariableSamples.col(0) *= 0.0;

    Eigen::MatrixXd posteriorSamples(stateSize, Nseeds);
    Eigen::MatrixXd priorSamples(stateSize, Nseeds);

    Eigen::VectorXd prior_mu = b * 0.0;
    Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(stateSize, stateSize) * 0.7;
    auto priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu, prior_cov); 
    for(int i = 0; i < Nseeds; i++)
    {
        priorSamples.col(i) = priorDensity->Sample();
    }
    posteriorSamples = priorSamples;
    
    sampleFlag = sampleDeltaStep;
    sampleI = 0;
    Eigen::MatrixXd forwardSamples(stateSize, Nseeds);
    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        std::cout << "Time " << time(timeI + 1) << std::endl;
        priorSamples = posteriorSamples;
        Eigen::MatrixXd forwardSamplesOld = forwardSamples;

        //Forecast step
        for(int i = 0; i < Nseeds; i++)
        {
            if(timeI == 0)
            {
	        forwardSamples.col(i) = A * x0 + priorDensity->Sample();//priorSamples.col(i);
            }
            else
            {
	        forwardSamples.col(i) = A * forwardSamplesOld.col(i) + priorDensity->Sample();//priorSamples.col(i);
            }
        }

        sampleFlag--;
        if(sampleFlag == 0)
        {
            sampleFlag = sampleDeltaStep;
	    Eigen::VectorXd meas = obs.col(sampleI);

	    //Kalman filter
	    posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(priorSamples, meas, meas_cov, H * forwardSamples);
            posteriorMean.col(timeI + 1) = posteriorSamples.rowwise().mean();

            for(int i = 0; i < Nseeds; i++)
            {
                if(timeI == 0)
                {
	            forwardSamples.col(i) = A * x0 + posteriorSamples.col(i);
                }
                else
                {
                    forwardSamples.col(i) = A * forwardSamplesOld.col(i)+ posteriorSamples.col(i);
                }
            }

            sampleI++;
        }
        else
        {
            for(int i = 0; i < Nseeds; i++)
            {
                posteriorSamples.col(i) += priorDensity->Sample();
            }
            posteriorMean.col(timeI + 1) = posteriorSamples.rowwise().mean();
        }

        singleVariableSamples.row(timeI + 1) = posteriorSamples.row(0); 
        stateRec.col(timeI + 1) = forwardSamples.rowwise().mean(); 
        minConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples, 0.05);
        maxConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples, 0.95);
    }

    ITHACAstream::exportMatrix(posteriorMean, "posteriorMean", "eigen", outputFolder);
    ITHACAstream::exportMatrix(obs, "obs", "eigen", outputFolder);
    ITHACAstream::exportMatrix(stateRec, "stateRec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(minConfidence, "minConfidence", "eigen", outputFolder);
    ITHACAstream::exportMatrix(maxConfidence, "maxConfidence", "eigen", outputFolder);
    ITHACAstream::exportMatrix(singleVariableSamples, "singleVariableSamples", "eigen", outputFolder);

    return 0;
}
