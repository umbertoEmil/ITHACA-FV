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
Class
    reducedInverseLaplacian
Description
    A reduced problem for the Laplace problem
SourceFiles
    reducedInverseLaplacian.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedInverseLaplacian class.
/// \dir
/// Directory containing the header and source files for the reducedInverseLaplacian class


#ifndef reducedInverseLaplacian_H
#define reducedInverseLaplacian_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "inverseLaplacianProblem.H"
#include <Eigen/Eigen>

/*---------------------------------------------------------------------------*\
                  Class reducedInverseLaplacian Declaration
\*---------------------------------------------------------------------------*/

/// Class to solve the online reduced problem associated with a Laplace problem
/** In this class are implemented the methods to solve a reduced version of the Laplace problem
    This class is used together with the laplacianProblem class */
class reducedInverseLaplacian: public reducedProblem
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedInverseLaplacian();

        /// Construct from an inverse Laplacian problem
        explicit reducedInverseLaplacian(inverseLaplacianProblem& problem);

        ~reducedInverseLaplacian() {};

        /// Problem object
        inverseLaplacianProblem* problem;

        /// Field of the measured temperatures derived from interpolation
        Eigen::MatrixXd Tmeas;

        /// Online computed reduced source terms
        Eigen::MatrixXd sourceT;
        Eigen::MatrixXd sourceLambda;
        Eigen::MatrixXd sourceDeltaT;

        /// Linear systems
        List<Eigen::MatrixXd> LinSysDirect;
        List<Eigen::MatrixXd> LinSysAdjoint;
        List<Eigen::MatrixXd> LinSysSensitivity;

        /// Reduced solutions
        Eigen::MatrixXd Treduced;
        Eigen::MatrixXd lambdaReduced;
        Eigen::VectorXd lambdaReducedSqr;
        Eigen::MatrixXd deltaTreduced;
        Eigen::MatrixXd TmeasReduced;

        /// Reconstructed fields
        List<Eigen::VectorXd> lambdaReducedList;
        List<scalar> gammaList;
        List<scalar> betaList;
        Eigen::VectorXd P; // Search direction
        Eigen::VectorXd Preduced;
        Eigen::VectorXd g; // Heat flux at GammaIn

        /// Offline computed matrices
        label offlineCG = 0;
        Eigen::MatrixXd L_T;
        Eigen::MatrixXd L_TinterpolationPlane;
        Eigen::MatrixXd L_lambda;
        Eigen::MatrixXd L_lambda_GammaIn;
        Eigen::MatrixXd L_deltaT;
        Eigen::MatrixXd F_T;
        Eigen::MatrixXd F_RBF;
        Eigen::MatrixXd B_deltaT;
        Eigen::MatrixXd b_deltaT_g;
        Eigen::MatrixXd faceArea_GammaIn;
        Eigen::VectorXd gammaOffline;
        Eigen::MatrixXd betaNum1;
        Eigen::MatrixXd betaNum2;
        Eigen::VectorXd betaDen;
        Eigen::MatrixXd Poffline;
        Eigen::MatrixXd PreducedOld;

        /// Conjugate Gradient parameters
        scalar J;
        scalar gamma;
        scalar beta;
        scalar betaOld;
        Eigen::VectorXd tempGamma;
        Eigen::MatrixXd tempGammaOld;

        /// Levenberg-Marquardt parameters
        Eigen::MatrixXd heatFluxBases;


        /// Online solution
        Eigen::MatrixXd online_solution;

        /// Counter for online sol
        label count_online_solve = 1;


        // Functions

        //--------------------------------------------------------------------------

        /// Set the BC and IF for the full sensitivity problem with homogeneous Neumann
        /// at the hotSide patch
        ///
        void setSensitivityBCandIF();

        //--------------------------------------------------------------------------

        /// Performs an online solution
        ///
        /// @param[in] field  Field to solve
        ///
        void assembleLinSys(volScalarField& field);

        //--------------------------------------------------------------------------

        /// Performs an online solution
        ///
        /// @param[in] field  Field to solve
        ///
        void solveOnline(volScalarField& field, label halfFull = 0);

        //--------------------------------------------------------------------------

        /// Fills the matrix Ldirect which contains the modes on the rows
        ///
        void get_L_T();

        //--------------------------------------------------------------------------

        /// Fills the matrix LdirectInterpolationPlane which contains the modes on the rows
        /// with zeros outside of the interpolation plane
        ///
        void get_L_TinterpolationPlane();

        //--------------------------------------------------------------------------

        /// Fills the matrix L_lambda which contains the modes on the rows
        ///
        void get_L_lambda();

        //--------------------------------------------------------------------------

        /// Fills the matrix L_deltaT which contains the modes on the rows
        ///
        void get_L_deltaT();

        //--------------------------------------------------------------------------

        ///  Assembles matrices F_T and TmeasReduced
        ///
        void assembleAdjointSourceReducedMatrices();

        //--------------------------------------------------------------------------

        ///  Assembles the offline component of the search direction Poffline
        ///
        void assembleG();

        //--------------------------------------------------------------------------

        ///  Assembles the matrices L_lambda_GammaIn, faceArea_GammaIn and B_deltaT
        ///
        void sensitivityBoundaryCoeff();

        //--------------------------------------------------------------------------

        ///  Assembles gammaOffline
        ///
        void conjugateCoeffOffline();

        //--------------------------------------------------------------------------

        ///  Assembles betaNum1, betaNum2 and betaDen
        ///
        void searchStepOffline();

        //--------------------------------------------------------------------------

        ///  Computes the conjugate coefficient, gamma
        ///
        void computeConjugateCoeff();

        //--------------------------------------------------------------------------

        ///  Compute the search steap, beta
        ///
        void searchStep();

        //--------------------------------------------------------------------------

        ///  Reconstruct the full order heat flux, g, and the search direction, P
        ///
        void reconstructHeatFlux();

        //--------------------------------------------------------------------------

        /// Reconstructs the full orded solution of the direct, adjoint and
        /// sensitivity problem. Saves the fields as member of problem.
        /// Converts g into a volScalarField, then it exports T, lambda, deltaT
        /// and g in the folder "folder" in a subdirectory named folderNumber+1
        ///
        /// @param[in] folderNumber  Fields are saved in the subdirectory folderNumber+1
        /// @param[in] folder        Name of the folder where the fields are saved
        ///
        void writeFields(label folderNumber, const char* folder);

        //--------------------------------------------------------------------------

        /// Computes the discrepancy functional, J, and check if convergence criteria are meet
        ///
        /// @param[out] 1 if one convergence criterium is meet, 2 if the algorithm
        ///             is diverging, 0 else
        ///
        int conjugateGradientConvergenceCheck();

        //--------------------------------------------------------------------------

        ///  Calls all the functions that asseble the matrices during the offline phase
        ///
        void conjugateGradientOffline();

        //--------------------------------------------------------------------------

        /// Reduced order conjugate gradient method
        ///
        /// @param[out]  return 1 if converged, 0 if reached maximum number of iterations
        ///
        int conjugateGradient();

        //--------------------------------------------------------------------------

        /// Conjugate gradient method in which only the linear systems are the reduced ones
        /// After the solution the full order fields are reconstructed
        ///
        /// @param[out]  return 1 if converged, 0 if reached maximum number of iterations
        ///
        int conjugateGradientHalfFull();

        //--------------------------------------------------------------------------

        ///
        void get_heatFluxBases();

        //--------------------------------------------------------------------------

        /// Reduced order Marquardt method
        ///
        /// @param[out]  return 1 if converged, 0 if reached maximum number of iterations
        ///
        int Marquardt();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






