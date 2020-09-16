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

    int Nseeds = 1000;

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

    std::cout << "We observe the state x by mean of the observation matrix \nH = " << H << std::endl; 
    std::cout << "The objective is to reconstruct the vector b kmowing H, x_0, and A" << std::endl; 
 

    int Ntimes = 201;
    int sampleDeltsStep = 3;
    double endTime = 10;
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(Ntimes, 0, endTime);

    Eigen::VectorXd xOld = x0;
    Eigen::MatrixXd X(stateSize,Ntimes);
    Eigen::MatrixXd forcer(stateSize,Ntimes);
    X.col(0) = x0;
    forcer.col(0) = b *0.0;

    int sampleFlag = sampleDeltsStep;
    int Nsamples = (Ntimes - 1) / sampleDeltsStep;
    int sampleI = 0;
    Eigen::MatrixXd obs(obsSize, Nsamples);

    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        forcer.col(timeI + 1) = - b * std::sin(time(timeI+1));
	Eigen::VectorXd xNew = A * xOld + forcer.col(timeI + 1);
	xOld = xNew;
	Eigen::VectorXd dNew = H * xNew;
        X.col(timeI + 1) = xNew;
        sampleFlag--;
        if(sampleFlag == 0)
        {
            sampleFlag = sampleDeltsStep;
            obs.col(sampleI) = dNew;
            sampleI++;
        }
    }
    std::cout << "\nobs = \n" << obs << std::endl;
    M_Assert(Nsamples == sampleI, "Somthing went wrong in the sampling");

    ITHACAstream::exportVector(time, "time", "eigen", outputFolder);
    ITHACAstream::exportMatrix(X, "X", "eigen", outputFolder);
    ITHACAstream::exportMatrix(forcer, "forcer", "eigen", outputFolder);

    Eigen::VectorXd x = x0;

    Eigen::MatrixXd prior_cov = Eigen::MatrixXd::Identity(stateSize, stateSize) * 0.1;


    Eigen::MatrixXd meas_cov = Eigen::MatrixXd::Identity(obsSize, obsSize) * 0.05;
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
    auto priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu, prior_cov); 
    for(int i = 0; i < Nseeds; i++)
    {
        priorSamples.col(i) = priorDensity->Sample();
    }
    posteriorSamples = priorSamples;
    
    sampleFlag = sampleDeltsStep;
    sampleI = 0;
    Eigen::VectorXd xNew = x0;
    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
        std::cout << "Time " << time(timeI + 1) << std::endl;
        xOld = xNew;
        sampleFlag--;
        if(sampleFlag == 0)
        {
            sampleFlag = sampleDeltsStep;
            priorSamples = posteriorSamples;
	    Eigen::VectorXd meas = obs.col(sampleI);

            //sampling
            Eigen::MatrixXd forwardSamples(stateSize, Nseeds);
            for(int i = 0; i < Nseeds; i++)
            {
                //priorSamples.col(i) += priorDensity->Sample();
	        forwardSamples.col(i) = A * xOld + priorSamples.col(i);
            }

	    //Kalman filter
	    posteriorSamples = ITHACAmuq::muq2ithaca::EnsembleKalmanFilter(priorSamples, meas, meas_cov, H * forwardSamples);

            singleVariableSamples.row(timeI + 1) = posteriorSamples.row(0); 

	    posteriorMean.col(timeI + 1) = posteriorSamples.rowwise().mean();

            sampleI++;
        }
        else
        {
            posteriorMean.col(timeI + 1) = posteriorMean.col(timeI);
        }

        xNew = A * xOld + posteriorMean.col(timeI + 1);
        stateRec.col(timeI + 1) = xNew; 
        minConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples, 0.05);
        maxConfidence.col(timeI + 1) = ITHACAmuq::muq2ithaca::quantile(posteriorSamples, 0.95);
    }


    ITHACAstream::exportMatrix(posteriorMean, "posteriorMean", "eigen", outputFolder);
    ITHACAstream::exportMatrix(stateRec, "stateRec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(minConfidence, "minConfidence", "eigen", outputFolder);
    ITHACAstream::exportMatrix(maxConfidence, "maxConfidence", "eigen", outputFolder);
    ITHACAstream::exportMatrix(singleVariableSamples, "singleVariableSamples", "eigen", outputFolder);

    return 0;
}
