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
    inverseLaplacianProblem_CG
Description
    A general full order implementation of an inverse problem
SourceFiles
    inverseLaplacianProblem_CG.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem_CG class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem_CG class.


#ifndef inverseLaplacianProblem_CG_H
#define inverseLaplacianProblem_CG_H
#include "inverseLaplacianProblem.H"
#define _USE_MATH_DEFINES

using namespace SPLINTER;

/// Class defining a plane for thermocouples measurements interpolation
class thermocouplesPlane
{
    public:
        double minX;
        double minZ;
        double Y;
        double maxX;
        double maxZ;
        List<scalar> Tmeas;
        List<scalar> Tdirect;
        List<scalar> Tsens;
        List<scalar> Tdiff;
        List<scalar> cellVol;
        List<int> cellID;
        List<scalar> thermocoupleX;
        List<scalar> thermocoupleZ;
        Foam::vector thermocoupleCellDim;
};

/// Class to implement a full order inverse parameterized laplacian problem
class inverseLaplacianProblem_CG: public inverseLaplacianProblem
{

    public:
        // Constructors
        inverseLaplacianProblem_CG();
        /// Construct with argc and argv
        inverseLaplacianProblem_CG(int argc, char* argv[]);
        virtual ~inverseLaplacianProblem_CG() {};


        /// Adjoint field
        autoPtr<volScalarField> _lambda;

        /// Sensibility temperature field
        autoPtr<volScalarField> _deltaT;

        /// Value fraction for the adjoint Robin BC
        List<scalar> valueFractionAdj;

        /// List of snapshots for the solutions
        PtrList<volScalarField> Tfield;
        PtrList<volScalarField> lambdaField;
        PtrList<volScalarField> deltaTfield;

        /// List of POD modes
        volScalarModes lambdaModes;
        volScalarModes deltaTmodes;

        /// Number of POD modes
        int NmodesT;
        int NmodesLambda;
        int NmodesDeltaT;


        int interpolation;
        int interpolationPlaneDefined;
        int cgIter;
        int cgIterMax;
        label saveSolInLists = 0;
        Eigen::MatrixXd Jlist; //cost funtion  [K^2]
        double gradJ_L2norm;
        double Jtol;
        double JtolRel;
        double gamma;
        double gamma_den;
        double beta;


        List<scalar> P;           //Search direction SHOULD BE [W/m2]
        List<scalar> gradJ;

        thermocouplesPlane interpolationPlane;
        Eigen::VectorXd Tsens;
        Eigen::MatrixXd muSamples;

        Eigen::VectorXd cellsInPlane;


        // Reduced Matrices
        /// A matrices
        List<Eigen::MatrixXd> ArbT;
        List<Eigen::MatrixXd> ArbLambda;
        List<Eigen::MatrixXd> ArbDeltaT;


        // Functions

        //--------------------------------------------------------------------------

        /// Set valueFraction list values for Robin condition
        void set_valueFraction();

        //--------------------------------------------------------------------------

        /// Set BC of the adjoint problem
        void assignAdjointBC();

        //--------------------------------------------------------------------------

        /// Assign the BC for the adjoint problem and returs the source field
        ///
        /// @param[out]  f  Source field of the adjoint problem
        ///
        volScalarField assignAdjointBCandSource();

        //--------------------------------------------------------------------------

        /// Set BC and IF of the sensitivity problem
        ///
        void assignSensitivityBC();

        //--------------------------------------------------------------------------

        /// Solve adjoint problem
        ///
        void solveAdjoint();

        //--------------------------------------------------------------------------

        /// Solve sensibility problem
        ///
        void solveSensitivity();

        //--------------------------------------------------------------------------

        /// Solve method for direct and sensitivity
        ///
        void solve(const char* problemID);

        //--------------------------------------------------------------------------

        /// Identifies the plane defined by the thermocouples
        ///
        void defineThermocouplesPlane();

        //--------------------------------------------------------------------------

        /// Fill the Foam::vector containing the values of the sensitivity solution at
        /// the thermocouples locations
        ///
        void sensibilitySolAtThermocouplesLocations();

        //--------------------------------------------------------------------------

        /// Conjugate gradient method
        ///
        /// @param[out]  1 if converged within cgIterMax iterations, 0 if not
        ///
        int conjugateGradient();

        //--------------------------------------------------------------------------

        /// Computes the gradient of cost function J and its L2 norm
        ///
        void computeGradJ();

        //--------------------------------------------------------------------------

        /// Computes the search direction P
        ///
        void searchDirection();

        //--------------------------------------------------------------------------

        /// Compute the search step beta
        ///
        void computeSearchStep();

        //--------------------------------------------------------------------------

        /// Updates the heat flux in the conjugate gradient iterations
        ///
        void updateHeatFlux();

        //--------------------------------------------------------------------------

        /// Convergence cher for the conjugate gradient method
        ///
        /// @param[out]  1 if convergence criteria met, 0 if not
        ///
        int conjugateGradientConvergenceCheck();

        //--------------------------------------------------------------------------

        /// Conjugate gradient method
        ///
        /// @param[in]  cx  X coordinate of the cell center
        /// @param[in]  cy  Y coordinate of the cell center
        /// @param[in]  cz  Z coordinate of the cell center
        /// @param[in]  thermocoupleCellDim  Dimension of the cell
        ///
        /// @param[out]  1 if converged within cgIterMax iterations, 0 if not
        ///
        int isInPlane(double cx, double cy, double cz,
                      Foam::vector thermocoupleCellDim);

        //--------------------------------------------------------------------------

        /// Writes fields to file
        ///
        /// @param[in]  folde  Folder in which the fields are written
        ///
        void writeFields(label folderNumber, const char* folder);

        //--------------------------------------------------------------------------

        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using radial basis functions
        ///
        void thermocouplesInterpolation();

        //--------------------------------------------------------------------------

        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using radial basis functions
        ///
        void thermocouplesInterpolation(DenseMatrix& RBFweights, DenseMatrix& RBFbasis);

        //--------------------------------------------------------------------------

        /// Computes the difference between direct problem solution and measure
        ///
        void differenceBetweenDirectAndMeasure();

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restart(word fieldName = "all");


        //--------------------------------------------------------------------------

        /// Performs offline solution
        ///
        //void offlineSolve();

};

#endif
