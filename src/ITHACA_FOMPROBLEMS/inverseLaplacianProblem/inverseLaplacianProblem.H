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
    inverseLaplacianProblem
Description
    A general full order implementation of an inverse problem
SourceFiles
    inverseLaplacianProblem.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the inverseLaplacianProblem class.
/// \dir
/// Directory containing the header and source files for the inverseLaplacianProblem class.


#ifndef inverseLaplacianProblem_H
#define inverseLaplacianProblem_H
#include <iostream>
#include "fvCFD.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
// #include "reducedLaplacian.H"
#include "ITHACAPOD.H"
#include "ITHACAutilities.H"
#include <Eigen/Dense>
#include <cmath>
#include "Foam2Eigen.H"
#include "cnpy.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "reductionProblem.H"
#include "Modes.H"
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
class inverseLaplacianProblem: public laplacianProblem
{

    public:
        // Constructors
        inverseLaplacianProblem();
        /// Construct with argc and argv
        inverseLaplacianProblem(int argc, char* argv[]);
        virtual ~inverseLaplacianProblem() {};

        ITHACAparameters* para;

        // Dummy variables to transform inverseLaplacianFoam into a class
        /// Temperature field
        autoPtr<volScalarField> _T;

        /// Adjoint field
        autoPtr<volScalarField> _lambda;

        /// Sensibility temperature field
        autoPtr<volScalarField> _deltaT;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// Time
        autoPtr<Time> _runTime;

        /// Dummy thermal conductivity with unitary value
        dimensionedScalar DT;

        /// Members for parameterized heat flux
	word folderOffline = "./ITHACAoutput/offlineParamBC/";
        PtrList<volScalarField> Tbasis;
        PtrList<volScalarField> Tad_base;
        
	Eigen::VectorXd residual;
        Eigen::MatrixXd Theta;
        Eigen::MatrixXd gPODmodes;
        Eigen::VectorXd addSol;


        /// List of snapshots for the solutions
        PtrList<volScalarField> Tfield;
        PtrList<volScalarField> lambdaField;
        PtrList<volScalarField> deltaTfield;

        /// List of POD modes
        Modes<scalar> lambdaModes;
        Modes<scalar> deltaTmodes;

        /// Number of POD modes
        int NmodesT;
        int NmodesLambda;
        int NmodesDeltaT;


        int interpolation;
        int interpolationPlaneDefined;
        int cgIter;
        int cgIterMax;
        int offlineSolutionI;
        int offlinePhase;
        int thermocouplesRead;
        int thermocouplesNum;
        label saveSolInLists = 0;
        double J; //cost funtion  [K^2]
        Eigen::MatrixXd Jlist; //cost funtion  [K^2]
        double gradJ_L2norm;
        double Jtol;
        double JtolRel;
        double k;  //thermal diffusivity [W/(m K)]
        double H;  //heat transfer coefficient [W/(m2 K)]
        double L2norm;
        double LinfNorm;
        double gamma;
        double gamma_den;
        double beta;
        double castWidth = 0;

        scalar v = 0.0;
        scalar homogeneousBC = 0.0;
        List<scalar> homogeneousBCcoldSide;
        List<scalar> Tf; //temperature at coldSide [K]
        List<scalar> refGrad;
        List<scalar> valueFraction;
        List<scalar> valueFractionAdj;

        label nProcs;
        label hotSide_ind;
        label coldSide_ind;

        List<scalar>
        g;           //Heat flux at hotSide for the direct problem [W/m2]
        List<List<scalar>> gBaseFunctions;
        List<scalar> gWeights;
        List<List<scalar>> gList;
        List<scalar> gTrue;       //True heat flux at hotSide [w/m2]
        List<scalar> faceCellArea; //[m2]
        List<scalar> P;           //Search direction SHOULD BE [W/m2]
        List<scalar> err;         //relative error []
        List<scalar> gradJ;


        // Marquardt method
        Eigen::MatrixXd Jacobian;
        Eigen::MatrixXd diagJacobian;
        label MarquardtHeatFluxSampling = 0;
        label MarquardtHeatFluxPOD = 0;

        thermocouplesPlane interpolationPlane;
        List<vector> thermocouplesPos;
        List<int> thermocouplesCellID;
        List<int> thermocouplesCellProc;
        List<Foam::vector> thermocouplesCellC;
        Eigen::VectorXd Tmeas;   //Temperature at the thermocouples locations [K]
        Eigen::VectorXd Tdirect;
        Eigen::VectorXd Tdiff;
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

        /// Set initial heat flux for the conjugate gradient method
        void set_g();

        //--------------------------------------------------------------------------

        /// Define the base functions used for the parametrization of g
        virtual void set_gBaseFunctions(word type, scalar shapeParameter = 1);

        //--------------------------------------------------------------------------

        /// Performs POD on the RBF
        void set_gBaseFunctionsPOD(label Nmodes);

        //--------------------------------------------------------------------------

        /// Set initial heat flux for the Levenberg-Marquardt method
        void set_gParametrized(word baseFuncType, scalar shapeParameter = 1);

        //--------------------------------------------------------------------------

        /// Convert list of boundary heat flux into field
        volScalarField list2Field(List<scalar> list, scalar innerField = 0.0);

        //--------------------------------------------------------------------------

        /// Convert Eigen matrix of boundary heat flux into field
        volScalarField eigen2Field(Eigen::MatrixXd matrix);

        //--------------------------------------------------------------------------

        /// Update the boundary condition g whe g is parameterized
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void update_gParametrized(List<scalar> weigths);

        //--------------------------------------------------------------------------

        ///
        ///
        Eigen::VectorXd TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b, label filter);

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
        void parameterizedBC(word linSys_solver = "fullPivLU",
                             label TSVD_filter = 0);

        //--------------------------------------------------------------------------

        ///
        ///
        void parameterizedBCpostProcess(Eigen::VectorXd weigths);

        //--------------------------------------------------------------------------

        /// Makes all the preliminary computations for the Marquardt method
        ///
        void MarquardtMethodSetUp();

        //--------------------------------------------------------------------------

        /// Update Jacobian matrix
        ///
        /// @param[in]  gWeightsOld    Weigths of parameterized BC at previous iteration
        /// @param[in]  weigthsUpdate  Update of the weigths in the present iteration
        /// @param[in]  Tcomp          Computed temperature at current iteration
        /// @param[in]  TcompOld       Computed temperature at previous iteration
        ///
        /// @param[out] Jacobian       Jacobian matrix
        ///
        Eigen::MatrixXd MarquardtUpdateJacobian(List<scalar> gWeightsOld,
                                                Eigen::VectorXd weigthsUpdate, Eigen::VectorXd Tcomp,
                                                Eigen::VectorXd TcompOld);


        //--------------------------------------------------------------------------

        /// Set valueFraction list values for Robin condition
        void set_valueFraction();

        //--------------------------------------------------------------------------

        /// Perform a solution which is used to fill the value of the temperature
        /// in the "virtual" thermocouples
        virtual void solveTrue();

        //--------------------------------------------------------------------------

        /// Set BC of the direct problem
        virtual void assignDirectBC();

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

        /// Set BC and IF of the additional problem for the parameterized heat flux
        ///
        virtual void solveAdditional();

        //--------------------------------------------------------------------------

        /// Solve direct problem
        ///
        void solveDirect();

        //--------------------------------------------------------------------------

        /// Solve adjoint problem
        ///
        void solveAdjoint();

        //--------------------------------------------------------------------------

        /// Solve sensibility problem
        ///
        void solveSensitivity();

        //--------------------------------------------------------------------------

        /// Solve laplacian problem without source term
        ///
        void solve(const char* problem);

        //--------------------------------------------------------------------------

        /// Identifies in the mesh the cells corresponding to the
        /// termocouples locations and creates the interpolation plane
        ///
        virtual void readThermocouples();

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[in]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(volScalarField& field);


        //--------------------------------------------------------------------------

        /// Identifies the plane defined by the thermocouples
        ///
        void defineThermocouplesPlane();

        //--------------------------------------------------------------------------

        /// Computes the difference between direct problem solution and measures
        ///
        void differenceBetweenDirectAndMeasure();

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

        /// Levenberg-Marquardt method
        ///
        /// @param[in]  exportSolutions 1 for saving the solutions at each iteration
        /// @param[in]  folder          Folder where the solutions are saved
        ///
        void MarquardtMethod(label exportSolutions, word folder,
                             label updateJacobian = 0);

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

        /// Compute maximum cell dimension in x, y and z
        ///
        /// @param[in]  ff
        /// @param[in]  pp
        /// @param[in]  cc
        /// @param[in]  pLabels
        /// @param[in]  pLocal
        ///
        /// @param[out] dim Maximum cell dimensions
        ///
        Foam::vector cellDim(const faceList& ff, const pointField& pp,
                             const cell& cc, labelList pLabels, pointField pLocal);

        //--------------------------------------------------------------------------

        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using radial basis functions
        ///
        void thermocouplesInterpolation(DenseMatrix& RBFweights, DenseMatrix& RBFbasis);

        //--------------------------------------------------------------------------

        /// Interpolates the thermocouples measurements in the plane defined
        /// in readThermocouples() using radial basis functions
        ///
        void thermocouplesInterpolation();

        //--------------------------------------------------------------------------

        /// Performs offline solution
        ///
        void offlineSolve();

        //--------------------------------------------------------------------------

        /// Samples the Tmeas and saves the heat fluxes
        ///
        void samplingHeatFluxMarquardt(word folder);

        //--------------------------------------------------------------------------

        /// Performs POD on the heat fluxes
        ///
        void heatFluxPodMarquardt(word folder);

        //--------------------------------------------------------------------------

        /// Offline sampling of direct problem with respect to heat flux parameters
        ///
        void MarquardtOffline(word folder);

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restart(word fieldName = "all");

        //--------------------------------------------------------------------------

        /// Reconsttuct online the field T 
        ///
        void reconstructT();

        //--------------------------------------------------------------------------

        ///// Offline sampling of direct problem with respect to heat flux parameters
        /////
        //void MarquardtOffline2(word folder);

};

#endif
