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
    EnKFtest.C 
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
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"


int main(int argc, char* argv[])
{
    word outputFolder = "./ITHACAoutput/";

    int stateSize = 4;
    int obsSize = 2;
    int Nseeds = 1000;

    Eigen::VectorXd x0(stateSize);
    x0 << 0.1, 0.2, 0.1, 0.2;

    Eigen::MatrixXd A(stateSize, stateSize);
    A << 1, 0, 0, 0,
	    0, 1.1, 0, 0,
	    0, 0, 1, .1,
	    0, 0, .1, 1;

    Eigen::MatrixXd H(obsSize, stateSize); //Observation Matrix
    H << 1, 0, 0, 0,
	    0, 0, 0, 1;
 
    Eigen::VectorXd b = x0 + Eigen::VectorXd::Ones(x0.size()) * 2.0;

    int Ntimes = 11;
    double endTime = 10;
    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(Ntimes, 0, endTime);

    Eigen::VectorXd xOld = x0;
    Eigen::MatrixXd obs(obsSize,Ntimes - 1);
    Eigen::MatrixXd X(stateSize,Ntimes);
    X.col(0) = x0;
    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
	Eigen::VectorXd xNew = A * xOld + b;
	xOld = xNew;
	Eigen::VectorXd dNew = H * xNew;
        X.col(timeI + 1) = xNew;
        obs.col(timeI) = dNew;
    }

    ITHACAstream::exportVector(time, "time", "eigen", outputFolder);
    ITHACAstream::exportMatrix(X, "X", "eigen", outputFolder);

    Eigen::VectorXd x = x0;

    Eigen::MatrixXd prior_cov(stateSize, stateSize);
    prior_cov << 0.5, .0, .0, .0,
           .0, 0.8, .0, .0,
	   .0, .0, .8, .0,
	   .0, .0, .0, .8;


    Eigen::MatrixXd meas_cov(obsSize, obsSize);
    meas_cov << 0.02, .0,
	   .0, .02;
    auto measNoise = std::make_shared<muq::Modeling::Gaussian>(Eigen::VectorXd::Zero(obsSize), meas_cov); 
    Eigen::MatrixXd posteriorMean(stateSize, Ntimes);
    posteriorMean.col(0) = x0;

    b = b *0.8;
    
    for(int timeI = 0; timeI < Ntimes - 1; timeI++)
    {
	Eigen::VectorXd meas = obs.col(timeI);
	Eigen::VectorXd prior_mu = posteriorMean.col(timeI);
        auto priorDensity = std::make_shared<muq::Modeling::Gaussian>(prior_mu, prior_cov); 

        //sampling
        Eigen::MatrixXd priorSamples(stateSize, Nseeds);
        Eigen::MatrixXd forwardSamples(stateSize, Nseeds);
        Eigen::MatrixXd measSamples(obsSize, Nseeds);
        for(int i = 0; i < Nseeds; i++)
        {
            priorSamples.col(i) = priorDensity->Sample();
	    forwardSamples.col(i) = A * priorSamples.col(i) + b;
	    measSamples.col(i) = meas + measNoise->Sample(); 
        }

	//Kalman filter
	Eigen::MatrixXd As = forwardSamples;
        for(int i = 0; i < Nseeds; i++)
        {
	    As.col(i) = forwardSamples.col(i) - forwardSamples.rowwise().mean();
        }

	Eigen::MatrixXd C = As * As.transpose() / (Nseeds - 1);

	Eigen::MatrixXd Y = measSamples - H * forwardSamples; //diff measurement data and simulated data

	Eigen::MatrixXd P = meas_cov + H * C * H.transpose(); //part which has to be inverted

	Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(P);
        auto P_rank = lu_decomp.rank();

	if (P_rank < P.cols() || P_rank < P.rows())
	{
	    Info << "P not invertible, exiting" << endl;
	    std::cout << P << std::endl;
	    exit(10);
	}
	else
	{
	    P = P.inverse();
	}

	Eigen::MatrixXd posteriorSample = forwardSamples + C * H.transpose() * P * Y;
	posteriorMean.col(timeI + 1) = posteriorSample.rowwise().mean();
    }
    ITHACAstream::exportMatrix(posteriorMean, "posteriorMean", "eigen", outputFolder);



     


    return 0;
}
