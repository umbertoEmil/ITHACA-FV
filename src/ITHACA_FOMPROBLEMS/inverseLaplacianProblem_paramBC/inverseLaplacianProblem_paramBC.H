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
    inverseLaplacianProblem_paramBC
Description
    A general full order implementation of an inverse problem
SourceFiles
    inverseLaplacianProblem_paramBC.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem_paramBC class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem_paramBC class.


#ifndef inverseLaplacianProblem_paramBC_H
#define inverseLaplacianProblem_paramBC_H

#include "inverseLaplacianProblem.H"

#define _USE_MATH_DEFINES

/// Class to implement a full order inverse laplacian problem
class inverseLaplacianProblem_paramBC: public inverseLaplacianProblem
{

    public:
        // Constructors
        inverseLaplacianProblem_paramBC();
        /// Construct with argc and argv
        inverseLaplacianProblem_paramBC(int argc, char* argv[]);
        virtual ~inverseLaplacianProblem_paramBC() {};


        /// Members for parameterized heat flux
	word folderOffline = "./ITHACAoutput/offlineParamBC/";
        PtrList<volScalarField> Tbasis;
        PtrList<volScalarField> Tad_base;
        
	Eigen::VectorXd residual;
        Eigen::MatrixXd Theta;
        Eigen::MatrixXd gPODmodes;
        Eigen::VectorXd addSol;

        int offlineSolutionI;
        int offlinePhase;

        List<List<scalar>> gBaseFunctions;
        List<scalar> gWeights;


        // Functions

        //--------------------------------------------------------------------------

        /// Set initial heat flux for the conjugate gradient method
        void set_g();

        //--------------------------------------------------------------------------

        /// Define the base functions used for the parametrization of g
        virtual void set_gBaseFunctions(word type, scalar shapeParameter = 1);

        //--------------------------------------------------------------------------

        /// Performs POD on the RBF
        void set_gBaseFunctionsPOD(label Nmodes);

        //--------------------------------------------------------------------------

        /// Set initial heat flux for the parameterized BC method
        void set_gParametrized(word baseFuncType, scalar shapeParameter = 1);

        //--------------------------------------------------------------------------

        /// Update the boundary condition g whe g is parameterized
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void update_gParametrized(List<scalar> weigths);

        //--------------------------------------------------------------------------

        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory ""./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        void parameterizedBCoffline(bool force = 0);

        //--------------------------------------------------------------------------

        ///
        ///
	Eigen::VectorXd parameterizedBC(word linSys_solver = "fullPivLU",
                             double regPar = 0);

        //--------------------------------------------------------------------------

        ///
        ///
        void parameterizedBCpostProcess(Eigen::VectorXd weigths);

        //--------------------------------------------------------------------------

        /// Set BC and IF of the additional problem for the parameterized heat flux
        ///
        virtual void solveAdditional();

        //--------------------------------------------------------------------------

        /// Reconsttuct online the field T 
        ///
        void reconstructT();
};

#endif
