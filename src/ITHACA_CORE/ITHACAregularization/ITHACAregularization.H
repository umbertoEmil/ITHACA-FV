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
    ITHACAregularization
Description
    Methods for regularizing ill-conditione linear systems
SourceFiles
    ITHACAregularization.C
    ITHACAregularizationTemplates.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAregularization class, it contains the implementation of
/// several methods for input output operations.
/// \dir
/// Directory containing the header, source and template files for the ITHACAregularization class, ITHACAparameters and Assert Functions.

#ifndef ITHACAregularization_H
#define ITHACAregularization_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdio>
#include <typeinfo>
#include <iostream>
#include <cassert>
#include <zlib.h>
#include <map>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <Eigen/Eigen>
#include "EigenFunctions.H"
#include "ITHACAutilities.H"
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"


#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop

#define MAXBUFSIZE (static_cast<int> (1e6))
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic pop

/// Namespace for regularization of ill-conditione linear system
namespace ITHACAregularization
{

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A      Matrix of coefficient
/// @param[in]  b      Column vector of solutions
/// @param[in]  filter Number of singular values to keep
///
/// @return     Column vector of variables
///
Eigen::VectorXd  TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b,
                      int filter);

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A                Matrix of coefficient
/// @param[in]  b                Column vector of solutions
/// @param[in]  parameterMethod  Regularization parameter selection method
///
/// @return     Column vector of variables
///
Eigen::VectorXd  TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b,
                      double noiseVariance, word parameterMethod);

//--------------------------------------------------------------------------
/// @brief      Truncated Singular Value regularization
///
/// @param[in]  A                Matrix of coefficient
/// @param[in]  b                Column vector of solutions
/// @param[in]  parameterMethod  Regularization parameter selection method
///
/// @return     Column vector of variables    
///
Eigen::VectorXd  TSVD(Eigen::MatrixXd A, Eigen::MatrixXd b, 
        double noiseVariance, word parameterMethod);

//--------------------------------------------------------------------------
/// @brief      Tikhonov regularization
///
/// @param[in]  A      Matrix of coefficient
/// @param[in]  b      Column vector of solutions
/// @param[in]  regularizationParameter Regularization parameter
///
/// @return     Column vector of variables
///
Eigen::VectorXd  Tikhonov(Eigen::MatrixXd A, Eigen::MatrixXd b,
                          double regularizationParameter);

//--------------------------------------------------------------------------
/// @brief      Tikhonov regularization
///
/// @param[in]  U      
/// @param[in]  s     
/// @param[in]  V      
/// @param[in]  b      Column vector of solutions
/// @param[in]  regularizationParameter Regularization parameter
///
/// @return     Column vector of variables
///
Eigen::VectorXd  Tikhonov(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd V, 
        Eigen::MatrixXd b, double regularizationParameter);

//--------------------------------------------------------------------------
/// @brief      Conjugate gradient applied implictly on the normal equation A'*A*x = A'*b
///
/// @param[in]  A      
/// @param[in]  b      Column vector of solutions
/// @param[in]  regularizationParameter Regularization parameter
///
/// @return     Column vector of variables
///
Eigen::VectorXd conjugateGradient(Eigen::MatrixXd A, Eigen::MatrixXd b, label Nsteps, bool reorth = 0, Eigen::VectorXd s = Eigen::VectorXd::Zero(1));

//--------------------------------------------------------------------------
/// @brief      Preconditioned conjugate gradient applied implictly on the normal equation A'*A*x = A'*b
///
/// @param[in]  A      
/// @param[in]  b      Column vector of solutions
/// @param[in]  regularizationParameter Regularization parameter
///
/// @return     Column vector of variables
///
Eigen::VectorXd PCGLS(Eigen::MatrixXd A, Eigen::MatrixXd L, Eigen::MatrixXd W, 
        Eigen::MatrixXd b, label Nsteps, bool reorth = 0, 
        Eigen::VectorXd s = Eigen::VectorXd::Zero(1));


//--------------------------------------------------------------------------
/// @brief      Generalized Cross Validation
///
/// @param[in]  U      
/// @param[in]  s     
/// @param[in]  b      Column vector of solutions
///
/// @return     Column vector of variables
///
scalar GCV(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b, scalar regPar0, word method);

struct GCVfun_properties
{
    Eigen::VectorXd s2;
    Eigen::VectorXd beta;
    scalar delta0;
    label MnN;
};
        
inline
double
GCVfun(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out,
          void* opt_data);


//--------------------------------------------------------------------------
/// @brief      L-curve
///
/// @param[in]  U      
/// @param[in]  s     
/// @param[in]  b      Column vector of solutions
///
/// @return     Column vector of variables
///
void Lcurve(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b, word method, label step = 0);

//--------------------------------------------------------------------------
/// @brief      L-curve
///
/// @param[in]  A      
/// @param[in]  b      Column vector of solutions
///
/// @return     Column vector of variables
///
void Lcurve_CG(Eigen::MatrixXd A, Eigen::MatrixXd b);

//--------------------------------------------------------------------------
/// @brief      Picard
///
/// @param[in]  U      
/// @param[in]  s     
/// @param[in]  b      Column vector of solutions
///
/// @return     Column vector of variables
///
void Picard(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b);


//--------------------------------------------------------------------------
/// @brief      Computes the discrete approximation L to the derivative operator
///             of order d on a regular grid with n points, i.e. L is (n-d)-by-n.
///
/// @param[in]  n       Grid size 
/// @param[in]  d       Derivative order
/// @param[in]  W       Orthonormal basis for the null space of L
///
/// @return     Column vector of variables
///
Eigen::MatrixXd get_L(label n, label d, Eigen::Ref<Eigen::MatrixXd> W);


//--------------------------------------------------------------------------
/// @brief      Method for calculating the pseudo-Inverse as recommended by Eigen developers 
///
/// @param[in]  
///
/// @return    
///
Eigen::MatrixXd pseudoInverse(Eigen::MatrixXd a, double epsilon = std::numeric_limits<double>::epsilon());

//--------------------------------------------------------------------------
/// @brief      Utility routine for "preconditioned" iterative methods.
///             Computes the vector
///                x = (L_p)'*y
///              where L_p is the A-weighted generalized inverse of L.
///             
///              Here, L is a p-by-n matrix, W holds a basis for the null space of L,
///              and T is a utility matrix which should be computed by routine pinit.
///             
///              Alternatively, L is square, and W and T are not needed.
///             
///              If W and T are not specified, then instead the routine computes
///                x = inv(L(:,1:p))'*y(1:p) .
///             
///              Notice that x and y may be matrices, in which case x(:,i)
///              corresponds to y(:,i).
///
/// @param[in]  
///
/// @return    
///
Eigen::VectorXd ltsolve(Eigen::MatrixXd L, Eigen::VectorXd y, 
        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(0,0), 
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(0,0));

//--------------------------------------------------------------------------
/// @brief      Utility routine for "preconditioned" iterative methods.
///             Computes the vector
///                x = (L_p)'*y
///              where L_p is the A-weighted generalized inverse of L.
///             
///              Here, L is a p-by-n matrix, W holds a basis for the null space of L,
///              and T is a utility matrix which should be computed by routine pinit.
///             
///              Alternatively, L is square, and W and T are not needed.
///             
///              If W and T are not specified, then instead the routine computes
///                x = inv(L(:,1:p))'*y(1:p) .
///             
///              Notice that x and y may be matrices, in which case x(:,i)
///              corresponds to y(:,i).
///
/// @param[in]  
///
/// @return    
///
Eigen::VectorXd lsolve(Eigen::MatrixXd L, Eigen::VectorXd y, 
        Eigen::MatrixXd W, Eigen::MatrixXd T);

//--------------------------------------------------------------------------
/// @brief      General method to precondition a linear system
///
/// @param[in]  linSys The linear system 
/// @param[in]  method Type of preconditioning to be used
///
/// @return    Preconditioned linear system
///
List<Eigen::MatrixXd> precondition(List<Eigen::MatrixXd> linSys, word method = "None", 
        int nSVD = 0);

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#endif
