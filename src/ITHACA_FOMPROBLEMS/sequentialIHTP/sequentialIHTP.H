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
    sequentialIHTP
Description
    A general full order implementation of an inverse problem
SourceFiles
    sequentialIHTP.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the sequentialIHTP class.
/// \dir
/// Directory containing the header and source files for the sequentialIHTP class.


#ifndef sequentialIHTP_H
#define sequentialIHTP_H
#include <iostream>
#include "fvCFD.H"
#include "interpolation.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IOmanip.H"
#include "Time.H"
#include "laplacianProblem.H"
// #include "reducedLaplacian.H"
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
class sequentialIHTP: public laplacianProblem
{

    public:
        // Constructors
        sequentialIHTP();
        /// Construct with argc and argv
        sequentialIHTP(int argc, char* argv[]);
        //sequentialIHTP(int argc, char* argv[], scalar _diffusivity);
        virtual ~sequentialIHTP() {};

        ITHACAparameters* para;

        // Dummy variables to transform inverseLaplacianFoam into a class
        /// Temperature field
        autoPtr<volScalarField> _T;

        /// List of temperature field in between two measurement instants
        PtrList<volScalarField> Ttime;

        /// List of temperature field in between two measurement instants
        PtrList<volScalarField> T_ic_time;

        /// Solutions for the space basis of the heat flux
        List<PtrList<volScalarField>> T_basis;

        /// Solutions of the direct problem with responsible of the initial condition
        List<PtrList<volScalarField>> T_IC;

        /// Mesh
        mutable autoPtr<fvMesh> _mesh;

        /// simpleControl
        autoPtr<simpleControl> _simple;

        /// fvOptions
        autoPtr<fv::options> _fvOptions;

        /// Time
        autoPtr<Time> _runTime;
        autoPtr<Time> _runTimeShort;

        /// Dummy thermal diffusivity with unitary value
        dimensionedScalar DT;

        /// Diffusivity value
	scalar diffusivity = 0.0; 

        /// Thermal conductivity [W/(m K)]
        double thermalCond;

        /// Heat transfer coefficient [W/(m2 K)]
        double HTC;

        /// Density [kg /m3]
        double density;

        /// Specific heat capacity [J/kg/K]
        double specificHeat;

        /// List of snapshots for the T_ic solutions
        PtrList<volScalarField> T_ic_field;

        /// List of POD modes
        volScalarModes T_ic_modes;

        /// T_ic reduced explicit matrix
        Eigen::MatrixXd T_ic_explicitMatrix_red;

        /// T_ic reduced implicit matrix
        Eigen::MatrixXd T_ic_implicitMatrix_red;

        /// Projection of the T_basis on the reduced space
        Eigen::MatrixXd T_basis_projectionMat;

        /// Magic points for the T_ic projection error estimation
        List<label> magicPoints;

        /// Matrix that reconstruct some points into the full order space
        Eigen::MatrixXd pointsReconstructMatrix;

        /// Matrix that reconstructs the T_basis at the magic points  
        Eigen::MatrixXd pointT_basis_reconstructionMat;

        /// L2 norm of the projection error of T_basis on the T_ic modes
        PtrList<volScalarField> projectionErrorT_basis;

        /// Number of POD modes
        int NmodesT_ic = 0;

        /// Time discretization (filled in the constructor)
        scalar startTime;
        scalar deltaTime;
        scalar endTime;
        label Ntimes;
        List<scalar> timeSteps;

        /// List of times at which the measurements are acquired
        /// (this List is filled by readThermocouples())
        List<scalar> samplingTime;

        /// Time interval in between samples (read from thermocouplesDict)
        scalar timeSamplesDeltaT;

        /// First sampling time (read from thermocouplesDict)
        scalar timeSamplesT0;

        /// Number of time samples (filled in createThermocouples.H)
        label timeSamplesNum;

        /// Number of timesteps between two samples
	label NtimeStepsBetweenSamples;

        /// Time sample index
	label timeSampleI;

        /// List of timesteps at which measurements are available
        List<label> samplingSteps; 

        /// 1 if readThermocouples() was called, 0 elsewise
        bool thermocouplesRead = 0;

        /// Number of thermocouples
        int thermocouplesNum;

        /// Cost funtion  [K^2]
        double J;

        /// List of cost funtions  [K^2]
        Eigen::VectorXd Jlist;

        /// Folder where the offline solutions are saved
        word folderOffline = "./ITHACAoutput/offlineParamBC_unsteady/";

        /// Folder where the offline solutions are saved
        word basisFolderOffline = folderOffline + "basis/";

        /// Number of timestep to solve for during offline phase
        label offlineTimestepsSize = 0;

        /// End time for the ofline computation
        scalar offlineEndTime = 0.0;

        bool offlineFlag = 0;

        /// Set to 1 when the Tb has been computed
        bool T_ic_ready = 0;

        /// Set to 1 for piecewise linear basis in time for the heat flux
        bool linearBasis = 0;

        scalar homogeneousBC = 0.0;

        /// Temperature at coldSide [K]
        List<scalar> Tf; 
        List<scalar> refGrad;
        List<scalar> valueFraction;

        /// Index of the hotSide patch
        label hotSide_ind;

        /// Index of the coldSide patch
        label coldSide_ind;

        /// Heat flux at hotSide, the first value is for t=0s
        List<List<scalar>> heatFlux;

        /// Heat flux space basis
        List<List<scalar>> heatFluxSpaceBasis;

        /// Time basis of the heat flux
        List<List<scalar>> heatFluxTimeBasis;

        /// Number of basis
        label Nbasis = 0;

        /// Number of samples considered in the offline phase
        label NsamplesWindow = 1;

        /// Weights of the parameterization
        List<scalar> heatFluxWeights;

        /// Weights of the parameterization
        List<scalar> heatFluxWeightsOld;

        /// Weights gradient of the parameterization
        List<scalar> heatFluxWeightsGrad;

        /// Weights gradient of the parameterization
        List<scalar> heatFluxWeightsGrad_old;

        /// True heat flux at hotSide [W/m2]
        List<List<scalar>> trueHeatFlux;

        /// Parametrized BC
        Eigen::VectorXd residual;

        /// Value of T_ic a the thermocouples
        Eigen::VectorXd T_ic_vector;
        
        /// Matrix for the inverse problem linear system
        Eigen::MatrixXd Theta;

        /// Type of solver used for the inverse problem linear system
        word linSys_solver;

        /// Number of singular values inc ase of TSVD regularization of the 
        /// inverse problem linear system
        label TSVD_filter;

        /// Parameter for the Tikhonov regularization of the inverse problem linear system
        scalar Tikhonov_filter;

        /// Number of steps for CG regularization of the inverse problem linear system
        label CG_Nsteps;


        /// Thermocouples positions
        List<vector> thermocouplesPos;
        
        /// ID of the cells containing the thermocouples
        List<int> thermocouplesCellID;

        /// Temperature at the thermocouples locations [K]
        Eigen::VectorXd Tmeas;   

        /// Temperature at the thermocouples locations in only two sampling times [K]
	Eigen::VectorXd TmeasShort;

        // Functions

        //--------------------------------------------------------------------------

        /// Set diffusivity
        void setDiffusivity(scalar _diff);

        //--------------------------------------------------------------------------

        /// Define the base functions used for the parametrization of g
        void setSpaceBasis(word type, scalar shapeParameter, label Npod = 0);

        //--------------------------------------------------------------------------

        /// Set parameterized heat flux defining the basis
        void setParametrizedHeatFlux(word spaceBaseFuncType,
            scalar shapeParameter_space, label Npod = 0);

        //--------------------------------------------------------------------------

        /// Update the boundary heat flux
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        void updateHeatFlux(List<scalar> weights);

        //--------------------------------------------------------------------------

        /// Update the boundary heat flux time basis
        ///
        /// @param[in]  weigths  New values of the weights of the basis functions
        ///
        virtual void updateTimeHeatFlux(List<scalar> weights) = 0;

        //--------------------------------------------------------------------------

        /// Performs offline computation for the parameterized BC method, if
        /// the offline directory ""./ITHACAoutputs/offlineParamBC" exists,
        /// it reads the solution from there
        ///
        /// @param[in]  force   If 1, force the offline phase to be computed
        ///
        virtual void offlinePhase(bool force = 0) = 0;

        //--------------------------------------------------------------------------

        /// Reconstructs the temperature field using superposition of effects
        ///
        virtual void reconstructT(word outputFolder) = 0;

        //--------------------------------------------------------------------------

        /// Reconstruct the field T at the input thermocouple 
        ///
        /// @param[in]   _TCindex  Index of the thermocouple
        ///
        /// @param[out]  Field value
        ///
        virtual scalar reconstructTatThermocouple(label _TCindex) = 0;

        //--------------------------------------------------------------------------

        /// Perform a solution of the direct problem
        ///
        /// @param[in]   _initialField  Initial field
        /// @param[in]   _outputFolder  Output folder
        ///
        void solveT(volScalarField _initialField, word _outputFolder);

        //--------------------------------------------------------------------------

        /// Compute the weights of the heat flux parameterization
        ///
        /// @param[in]   _outputFolder      Output folder
        /// @param[in]   _initialField      Initial field
        /// @param[in]   _initialWeights    Initial weights
        ///
        virtual void computeHeatFluxWeights(volScalarField _initialField, 
                List<scalar> _initialWeights, word _outputFolder) = 0;

        //--------------------------------------------------------------------------

        /// Perform a solution of a linear system
        /// @param[in]   _linSys    Linear system matrix and source vector
        ///
        /// @param[out]  Solution vector
        ///
        Eigen::VectorXd solveLinSys(List<Eigen::MatrixXd> _linSys); 

        //--------------------------------------------------------------------------

        /// Set valueFraction list values for Robin condition
        void set_valueFraction();

        //--------------------------------------------------------------------------

        /// Set BC of the direct problem
        virtual void assignDirectBC(label timeI);

        //--------------------------------------------------------------------------

        /// Solve the T_ic problem
        ///
        void solveT_ic(volScalarField _initialField);

        //--------------------------------------------------------------------------

        /// Compute T_ic modes prome snapshots
        ///
        void getT_ic_modes();

        //--------------------------------------------------------------------------

        /// Project T_ic matrices onto the reduced spaced
        ///
        void projectT_ic();

        //--------------------------------------------------------------------------

        /// Assemble the matrices to go from the gWeights to the T_ic_ reduced space
        ///
        void projectDirectOntoT_ic();

        //--------------------------------------------------------------------------

        /// Assemble the matrix pointsProjectionMatrix to project some points on the 
        /// reduced basis space
        ///
        void pointProjectionOffline(); 

        //--------------------------------------------------------------------------

        /// Compute the L2 norm of the projection error for each T_basis
        ///
        void projectionErrorOffline();

        //--------------------------------------------------------------------------

        /// Assemble all the matrices required in the online phase
        ///
        void T_ic_offline(int NmagicPoints);

        //--------------------------------------------------------------------------

        /// Solve direct problem
        ///
        void solveDirect();

        //--------------------------------------------------------------------------

        /// Identifies in the mesh the cells corresponding to the
        /// termocouples locations.
        ///
        virtual void readThermocouples();

        //--------------------------------------------------------------------------

        /// Interpolates the field value at a point
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[out]  Vector of field values at thermocouples points
        ///
        scalar fieldValueAtPoint(volScalarField& _field, vector _point);

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        /// NOTE: do NOT call when field is an element of a PtrList
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[out]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(volScalarField field);

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        /// Always use thisone for elementes of PtrList
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[out]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(PtrList<volScalarField> fieldList,
                label fieldI);

        //--------------------------------------------------------------------------

        /// Interpolates the field value at the thermocouples points
        ///
        /// @param[in]  field    Field to read the values
        ///
        /// @param[out]  Vector of field values at thermocouples points
        ///
        Eigen::VectorXd fieldValueAtThermocouples(PtrList<volScalarField> fieldList);

        //--------------------------------------------------------------------------

        /// Restart temperature field
        ///
        void restart();

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restartOffline();

        //--------------------------------------------------------------------------

        /// Restart fields
        ///
        void restartT_ic();

        //--------------------------------------------------------------------------

        /// Fills the vector samplingSteps which contains the timesteps at which the 
        /// measurements are taken
        ///
        void sampling2symulationTime();

        //--------------------------------------------------------------------------

        ///
        ///
        void parameterizedHeatFlux_postProcess(List<Eigen::MatrixXd> linSys,
                                         Eigen::VectorXd weigths, 
                                         word outputFolder, label verbose = 0);

        //--------------------------------------------------------------------------

        /// Find the points at with the projection error is computed
        ///
        void findMagicPoints(int NmagicPoints);

        //--------------------------------------------------------------------------

        /// Convert list of boundary heat flux into field
        volScalarField list2Field(List<scalar> list, scalar innerField = 0.0);

};

#endif
