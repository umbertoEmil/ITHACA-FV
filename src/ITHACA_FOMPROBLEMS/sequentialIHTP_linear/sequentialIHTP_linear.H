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
    sequentialIHTP_linear
Description
    A general full order implementation of an inverse problem
SourceFiles
    sequentialIHTP_linear.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the sequentialIHTP_linear class.
/// \dir
/// Directory containing the header and source files for the sequentialIHTP_linear class.


#ifndef sequentialIHTP_linear_H
#define sequentialIHTP_linear_H
#include <iostream>
#include "fvCFD.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "sequentialIHTP.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <cmath>
#include "Foam2Eigen.H"
#include "cnpy.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "ITHACAregularization.H"

#define _USE_MATH_DEFINES

using namespace SPLINTER;

/// Class to implement a full order inverse parameterized laplacian problem
class sequentialIHTP_linear: public sequentialIHTP
{

    public:
        // Constructors
        sequentialIHTP_linear();
        /// Construct with argc and argv
        sequentialIHTP_linear(int argc, char* argv[]);
        //sequentialIHTP_linear(int argc, char* argv[], scalar _diffusivity);
        virtual ~sequentialIHTP_linear() {};

        /// Solutions of the T_d problem
        List<PtrList<volScalarField>> T_d;

        /// Matrix for the inverse problem linear system
        Eigen::MatrixXd Theta_d;

        scalar onlineCPUtime_avg;

        word reconstructedTemperatureName = "recTemperature";
        word reconstructedHeatFluxName = "recHeatFluc";

        // Functions

        //--------------------------------------------------------------------------

        /// Update the boundary heat flux time basis
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void updateTimeHeatFlux(List<scalar> weights);

        //--------------------------------------------------------------------------

        /// Solve the T_d problem
        ///
        void solveT_d(label _baseI);

        //--------------------------------------------------------------------------

        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory ""./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        void offlinePhase(bool force = 0);

        //--------------------------------------------------------------------------

        /// Reconstructs the temperature field using superposition of effects
        ///
        void reconstructT(word _outputFolder);


        //--------------------------------------------------------------------------

        /// 
        ///
        scalar reconstructTatThermocouple(label _TCindex); 

        //--------------------------------------------------------------------------

        /// Compute the weights of the heat flux parameterization
        ///
        /// @param[in]   _outputFolder      Output folder
        /// @param[in]   _initialField      Initial field
        /// @param[in]   _initialWeights    Initial weights
        ///
        void computeHeatFluxWeights(volScalarField _initialField, 
                List<scalar> _initialWeights, word _outputFolder);
};

#endif
