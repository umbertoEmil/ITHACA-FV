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
#include "ITHACAregularization.H"
#include <Eigen/Dense>
#include <cmath>
#include "Foam2Eigen.H"
#include "cnpy.H"
#include "mixedFvPatchFields.H"
#include "cellDistFuncs.H"
#include "reductionProblem.H"
#include "Modes.H"
#define _USE_MATH_DEFINES

/// Class to implement a full order inverse laplacian problem
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

        double k;  //thermal diffusivity [W/(m K)]
        double H;  //heat transfer coefficient [W/(m2 K)]

        int thermocouplesRead;
        int thermocouplesNum;
        double J; //cost funtion  [K^2]
        double L2norm;
        double LinfNorm;

        /// Lists for the Robin BC
        scalar homogeneousBC = 0.0;
        List<scalar> homogeneousBCcoldSide;
        List<scalar> Tf; //temperature at coldSide [K]
        List<scalar> refGrad;
        List<scalar> valueFraction;

        label nProcs;
        label hotSide_ind;
        label coldSide_ind;

        List<scalar>
        g;           //Heat flux at hotSide for the direct problem [W/m2]
        List<List<scalar>> gList;
        List<scalar> gTrue;       //True heat flux at hotSide [w/m2]
        List<scalar> faceCellArea; //[m2]
        List<scalar> err;         //relative error []

        List<vector> thermocouplesPos;
        List<int> thermocouplesCellID;
        List<int> thermocouplesCellProc;
        List<Foam::vector> thermocouplesCellC;
        Eigen::VectorXd Tmeas;   //Temperature at the thermocouples locations [K]
        Eigen::VectorXd Tdirect;
        Eigen::VectorXd Tdiff;



        // Functions

        //--------------------------------------------------------------------------
        
        /// Set initial heat flux for the conjugate gradient method
        void set_g();

        //--------------------------------------------------------------------------

        /// Convert list of boundary heat flux into field
        volScalarField list2Field(List<scalar> list, scalar innerField = 0.0);

        //--------------------------------------------------------------------------

        /// Convert Eigen matrix of boundary heat flux into field
        volScalarField eigen2Field(Eigen::MatrixXd matrix);

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

        /// Solve direct problem
        ///
        void solveDirect();

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

        /// Computes the difference between direct problem solution and measures
        /// Saves the difference vector in Tdiff
        ///
        void differenceBetweenDirectAndMeasure();

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

        /// Restart fields
        ///
        void restart(word fieldName = "all");
};

#endif
