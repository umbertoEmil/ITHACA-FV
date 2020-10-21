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
    reducedSequentialIHTP
Description
    A reduced problem for the stationary NS equations
SourceFiles
    reducedSequentialIHTP.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the reducedSequentialIHTP class.
/// \dir
/// Directory containing the header and source files for the reducedSequentialIHTP class

#ifndef ReducedSequentialIHTP_H
#define ReducedSequentialIHTP_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedProblem.H"
#include "sequentialIHTP.H"
#include "ITHACAutilities.H"
#include "EigenFunctions.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>


/*---------------------------------------------------------------------------*\
                        Class reducedSequentialIHTP Declaration
\*---------------------------------------------------------------------------*/

/// Class where it is implemented a reduced problem for the <b> sequential inverse heat transfer </b> problem.
/** In this class are implemented the methods to solve a reduced version of the sequential inverse heat transfer problem 
equations. This class must be used together with the sequentialIHTP class */
class reducedSequentialIHTP: public sequentialIHTP
{
    private:

    public:
        // Constructors
        /// Construct Null
        reducedSequentialIHTP();
        reducedSequentialIHTP(int argc, char* argv[]);

        ~reducedSequentialIHTP() {};

        Eigen::VectorXd T0red;
        double T0projectionTol = 0;
        bool previousWasReduced = 0;
        Eigen::VectorXi onlineCountVec; 
        Eigen::VectorXd onlineWindowsVec; 


        // Functions
        
        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void parameterizedBC(word outputFolder, volScalarField initialField,
                Eigen::VectorXi errorCells, word linSys_solver = "fullPivLU", 
                label TSVD_filter = 0);

        /// Method to perform an online solve using a PPE stabilisation method
        ///
        /// @param[in]  vel_now  The vector of online velocity. It is defined in
        /// with an Eigen::MatrixXd and must have one col and as many rows as the number
        /// of parametrized boundary conditions.
        ///
        void solveT0online(volScalarField initialField, bool useReducedInitialField = 0);

        ///
        /// 
        ///
        /// @param[in]  
        /// @param[in]  
        ///
        double T0projectionError(volScalarField& T0in, Eigen::VectorXi errorCells);

        ///
        /// 
        ///
        volScalarField reconstrucT_lastTime();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif






