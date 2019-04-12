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
     foam2eigen
 Description
     Converts foam objects into Eigen and viceversa
 SourceFiles
     foam2eigen.C
 \*---------------------------------------------------------------------------*/

/// \file
/// Header file of the Foam2Eigen class.
/// \dir
/// Directory containing the header and source files for the Foam2Eigen class.

#ifndef Foam2Eigen_H
#define Foam2Eigen_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ITHACAassert.H"
#include <tuple>
#include <sys/stat.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
  Class foam2eigen Declaration
\*---------------------------------------------------------------------------*/

/// Class to converts OpenFOAM objects into Eigen and viceversa
class Foam2Eigen
{
    private:

    public:

        //----------------------------------------------------------------------
        /// @brief      Convert a FvMatrix OpenFOAM matrix (Linear System) into
        ///             a Eigen Matrix A and a source vector b
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] A                 The matrix or the list of matrices,
        ///                               can be dense or sparse
        /// @param[out] b                 The source term vector or the list of
        ///                               source term vectors, always dense
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be scalar
        ///                               or vector
        /// @tparam     type_A            The type of matrix or the list of
        ///                               matrices, can be dense
        ///                               (Eigen::MatrixXd) or sparse
        ///                               (Eigen::SparseMatrix<double>) or a
        ///                               List of them
        /// @tparam     type_B            The type source term vector or the
        ///                               list can be Eigen::VectorXd or
        ///                               List<Eigen::VectorXd>
        ///
        template <class type_foam_matrix, class type_A, class type_B>
        static void fvMatrix2Eigen(fvMatrix<type_foam_matrix> foam_matrix, type_A& A,
                                   type_B& b);

        //----------------------------------------------------------------------
        /// @brief      Convert a ldu OpenFOAM matrix into a Eigen Matrix A
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] A                 The matrix or the list of matrices,
        ///                               can be dense or sparse
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be
        ///                               fvScalarMatrix or fvVectorMatrix
        /// @tparam     type_A            The type of matrix or the list of
        ///                               matrices, can be dense
        ///                               (Eigen::MatrixXd) or sparse
        ///                               (Eigen::SparseMatrix<double>) or a
        ///                               List of them
        ///
        template <class type_foam_matrix, class type_A>
        static void fvMatrix2EigenM(fvMatrix<type_foam_matrix>& foam_matrix, type_A& A);

        //----------------------------------------------------------------------
        /// @brief      Convert a ldu OpenFOAM matrix into a source vector b
        ///
        /// @param[in]  foam_matrix       The foam matrix can be fvScalarMatrix
        ///                               or fvVectorMatrix
        /// @param[out] b                 The source term vector or the list of
        ///                               source term vectors, always dense
        ///
        /// @tparam     type_foam_matrix  The type of foam matrix can be
        ///                               fvScalarMatrix or fvVectorMatrix
        /// @tparam     type_B            The type source term vector or the
        ///                               list can be Eigen::VectorXd or
        ///                               List<Eigen::VectorXd>
        ///
        template <class type_foam_matrix, class type_B>
        static void fvMatrix2EigenV(fvMatrix<type_foam_matrix>& foam_matrix, type_B& b);


        //----------------------------------------------------------------------
        /// @brief      Convert a PtrList of snapshots to Eigen matrix (only
        ///             internal field)
        ///
        /// @param[in]  fields      The fields can be a PtrList<volScalarField>
        ///                         or PtrList<volVectorField>
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     An Eigen matrix containing as columns the snapshots
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::MatrixXd PtrList2Eigen(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief      Convert a vector OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd field2Eigen(GeometricField<vector, PatchField, GeoMesh>&
                                           field);

        //----------------------------------------------------------------------
        /// @brief      Convert a scalar OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd field2Eigen(GeometricField<scalar, PatchField, GeoMesh>&
                                           field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field    The field
        ///
        /// @tparam     Type     vector or scalar.
        /// @tparam     GeoMesh  volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<class Type, class GeoMesh>
        static Eigen::VectorXd field2Eigen(const
                                           DimensionedField<Type, GeoMesh>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM field into an Eigen Vector
        ///
        /// @param[in]  field  The field
        ///
        /// @tparam     Type   Type of the field can be volScalarField or
        ///                    volVectorField
        ///
        /// @return     Dense   Eigen vector
        ///
        template<class Type>
        static Eigen::VectorXd field2Eigen(const Field<Type>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM scalar field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::VectorXd> field2EigenBC(
            GeometricField<scalar, PatchField, GeoMesh>& field);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM vector field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  field       The field
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen vector
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::VectorXd> field2EigenBC(
            GeometricField<vector, PatchField, GeoMesh>& field);


        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM scalar field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  fields      The field
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen matrix
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::MatrixXd> PtrList2EigenBC(
            PtrList<GeometricField<scalar, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief      Convert an OpenFOAM vector field to a List of Eigen
        ///             Vectors, one for each boundary
        ///
        /// @param[in]  fields      The field
        /// @param[in]  Nfields     The nfields
        ///
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Dense Eigen matrix
        ///
        template<template<class> class PatchField, class GeoMesh>
        static List<Eigen::MatrixXd> PtrList2EigenBC(
            PtrList<GeometricField<vector, PatchField, GeoMesh>>& fields,
            label Nfields = -1);

        //----------------------------------------------------------------------
        /// @brief         Convert a vector in Eigen format into an OpenFOAM
        ///                scalar GeometricField
        ///
        /// @param[in/out] field         OpenFOAM GeometricField
        /// @param[in]     eigen_vector  Vector in Eigen format
        ///
        /// @tparam        PatchField    fvPatchField or fvsPatchField.
        /// @tparam        GeoMesh       volMesh or surfaceMesh.
        ///
        /// @return        OpenFOAM GeometricField
        ///
        template<template<class> class PatchField, class GeoMesh>
        static GeometricField<scalar, PatchField, GeoMesh> Eigen2field(
            GeometricField<scalar, PatchField, GeoMesh>& field,
            Eigen::VectorXd& eigen_vector);

        //----------------------------------------------------------------------
        /// @brief         Convert a vector in Eigen format into an OpenFOAM
        ///                vector GeometricField
        ///
        /// @param[in/out] field         OpenFOAM GeometricField
        /// @param[in]     eigen_vector  Vector in Eigen format
        ///
        /// @tparam        PatchField    fvPatchField or fvsPatchField.
        /// @tparam        GeoMesh       volMesh or surfaceMesh.
        ///
        /// @return        OpenFOAM GeometricField
        ///
        template<template<class> class PatchField, class GeoMesh>
        static GeometricField<vector, PatchField, GeoMesh> Eigen2field(
            GeometricField<vector, PatchField, GeoMesh>& field,
            Eigen::VectorXd& eigen_vector);

        //----------------------------------------------------------------------
        /// @brief         Converts a matrix in Eigen format into an OpenFOAM
        ///                Field
        ///
        /// @param[in/out] field   OpenFOAM Field
        /// @param[in]     matrix  Matrix in Eigen format
        ///
        /// @tparam        Type    vector or scalar.
        ///
        /// @return        OpenFOAM Field
        ///
        template <class Type>
        static Field<Type> Eigen2field(Field<Type>& field, Eigen::MatrixXd& matrix);

        //----------------------------------------------------------------------
        /// @brief      Obtain the Mass Matrix from a list of reduced basis
        ///
        /// @param      matrix      The matrix
        /// @param[in]  modes       The OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      Number of requested modes
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     An Eigen matrix containing the mass matrix
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static std::tuple<Eigen::MatrixXd, Eigen::VectorXd> projectFvMatrix(
            fvMatrix<Type>& matrix,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes);


        //----------------------------------------------------------------------
        /// @brief      Perform the projection of an OpenFOAM field onto a set
        ///             of modes using the Eigen matrix multiplication
        ///
        ///             It is the standard Galerkin projection of a field onto a
        ///             reduced basis:
        ///
        ///             @f[  \mathbf{f} \in R^{Nx1} , \mathbf{f_r} \in R^{nx1} ,
        ///             \mathbf{M} \in R^{Nxn}
        /// @f]
        ///
        ///             @f[ \mathbf{f_r} = \mathbf{M^T f}
        /// @f]
        ///
        /// @param[in]  field       The f OpenFOAM field
        /// @param[in]  modes       The M OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      The n number of modes you want to use (The
        ///                         size of input modes can be bigger than the
        ///                         requested number of modes)
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Eigen vector with the projected field
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::VectorXd projectField(
            GeometricField<Type, PatchField, GeoMesh>& field,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
            label Nmodes);

        //----------------------------------------------------------------------
        /// @brief      Obtain the Mass Matrix from a list of reduced basis
        ///
        /// @param      modes       The OpenFOAM modes contained in a PtrList
        /// @param[in]  Nmodes      The number of requested modes
        ///
        /// @tparam     Type        vector or scalar.
        /// @tparam     PatchField  fvPatchField or fvsPatchField.
        /// @tparam     GeoMesh     volMesh or surfaceMesh.
        ///
        /// @return     Eigen matrix containing the mass matrix
        ///
        template<class Type, template<class> class PatchField, class GeoMesh>
        static Eigen::MatrixXd MassMatrix(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes, label Nmodes);

        //----------------------------------------------------------------------
        /// @brief      Convert a PtrList of OpenFOAM fvMatrix into a tuple of
        ///             lists of Eigen Sparse Matrices and source vectors
        ///
        /// @param[in]  MatrixList  The matrix list
        ///
        /// @tparam     Type        { description }
        /// @tparam     [in]  Type  type of fvMatrix can be either scalar or vectorial
        ///
        /// @return     a tuple containing a list of sparse matrices and a list
        ///             of source vectors
        ///
        template<class Type>
        static std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
                LFvMatrix2LSM(PtrList<fvMatrix<Type>>& MatrixList);

        //--------------------------------------------------------------------------
        /// @brief      Convert a Foam List into an Eigen matrix with one column
        ///
        /// @param[in]  list        The Foam List
        ///
        /// @tparam     type_matrix The type of data in the List (int or double)
        ///
        /// @return     The Eigen matrix
        ///
        template <class type_matrix>
        static Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic>
        List2EigenMatrix ( List<type_matrix> list );

        //--------------------------------------------------------------------------
        /// @brief      Convert an Eigen matrix with one column into a Foam List
        ///
        /// @param[in]  matrix       The Eigen matrix
        ///
        /// @tparam     type_matrix The type of data in the List (int or double)
        ///
        /// @return     The Foam List
        ///
        template <class type_matrix>
        static List<type_matrix> EigenMatrix2List (
            Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix );



};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class type_f>
Eigen::VectorXd Foam2Eigen::field2Eigen(const
                                        DimensionedField<type_f, Foam::volMesh>& field)
{
    Eigen::VectorXd out;
    out.resize(int(field.size()));

    for (int l = 0; l < field.size(); l++)
    {
        out(l) = field[l];
    }

    return out;
}

template <>
Field<scalar> Foam2Eigen::Eigen2field(Field<scalar>& field,
                                      Eigen::MatrixXd& matrix);

template <>
Field<vector> Foam2Eigen::Eigen2field(Field<vector>& field,
                                      Eigen::MatrixXd& matrix);

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<vector, fvPatchField, volMesh>& field);

template<>
Eigen::VectorXd Foam2Eigen::field2Eigen(
    GeometricField<scalar, fvPatchField, volMesh>& field);

template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<vector, fvPatchField, volMesh>& field);

template<>
List<Eigen::VectorXd> Foam2Eigen::field2EigenBC(
    GeometricField<scalar, fvPatchField, volMesh>& field);

template<class type_f>
List<Eigen::MatrixXd> PtrList2EigenBC(type_f& fields, int Nfields);

template<>
List<Eigen::MatrixXd> PtrList2EigenBC(
    PtrList<GeometricField<scalar, fvPatchField, volMesh>>& fields,
    int Nfields);

template<>
List<Eigen::MatrixXd> PtrList2EigenBC(
    PtrList<GeometricField<vector, fvPatchField, volMesh>>& fields,
    int Nfields);

template<class type_matrix>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        Foam2Eigen::LFvMatrix2LSM(PtrList<fvMatrix<type_matrix>>& MatrixList)
{
    List<Eigen::SparseMatrix<double>> SM_list;
    List<Eigen::VectorXd> V_list;
    int LSize =  MatrixList.size();
    SM_list.resize(LSize);
    V_list.resize(LSize);
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;

    for (int j = 0; j < LSize; j++)
    {
        fvMatrix2Eigen(MatrixList[j], A, b);
        SM_list[j] = A;
        V_list[j] = b;
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(SM_list, V_list);
    return tupla;
}

template<class type_matrix>
Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic>
Foam2Eigen::List2EigenMatrix ( List<type_matrix> list )
{
    Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix(list.size(),
            1);

    for (label i = 0; i < matrix.rows(); i++)
    {
        matrix(i, 0) = list[i];
    }

    return matrix;
}

template<class type_matrix>
List<type_matrix> Foam2Eigen::EigenMatrix2List (
    Eigen::Matrix<type_matrix, Eigen::Dynamic, Eigen::Dynamic> matrix )
{
    if (matrix.cols() == 1)
    {
        List<int> list(matrix.rows());

        for (label i = 0; i < matrix.rows(); i++)
        {
            list[i] = matrix(i, 0);
        }

        return list;
    }
    else
    {
        Info << "I still must code EigenMatrix2List for matrix cols>1, exiting" << endl;
        exit(11);
    }
}

template <class type_m>
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> Foam2Eigen::projectFvMatrix(
    fvMatrix<type_m>& matrix,
    PtrList<GeometricField<type_m, fvPatchField, volMesh>>& modes, int Nmodes)
{
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd Ar;
    Eigen::VectorXd b;
    Eigen::VectorXd br;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Foam2Eigen::fvMatrix2Eigen(matrix, A, b);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }

    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
        }
    }

    Ar = Eig_Modes.transpose() * A * Eig_Modes;
    br = Eig_Modes.transpose() * b;
    std::tuple <Eigen::MatrixXd, Eigen::VectorXd> tupla;
    tupla = std::make_tuple(Ar, br);
    return tupla;
}

template <class type_PtrList>
Eigen::MatrixXd Foam2Eigen::MassMatrix(
    PtrList<GeometricField<type_PtrList, fvPatchField, volMesh>>& modes, int Nmodes)
{
    Eigen::MatrixXd Mr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), Nmodes);

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), Nmodes);
    }

    if (Volumes.rows() == Eig_Modes.rows())
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i) = Volumes;
        }
    }
    else
    {
        for (auto i = 0; i < Nmodes; i++)
        {
            VolumesN.col(i).segment(0, Volumes.rows()) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
            VolumesN.col(i).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
        }
    }

    Mr = Eig_Modes.transpose() * (Eig_Modes.cwiseProduct(VolumesN));
    return Mr;
}

template <class type_f>
Eigen::VectorXd Foam2Eigen::projectField(
    GeometricField<type_f, fvPatchField, volMesh>& field,
    PtrList<GeometricField<type_f, fvPatchField, volMesh>>& modes,
    int Nmodes)
{
    Eigen::VectorXd fr;
    Eigen::MatrixXd Eig_Modes = PtrList2Eigen(modes, Nmodes);
    Eigen::VectorXd f = Foam2Eigen::field2Eigen(field);
    Eigen::VectorXd Volumes = field2Eigen(modes[0].mesh());
    Eigen::MatrixXd VolumesN(Volumes.rows(), 1);
    VolumesN = Volumes;

    if (Volumes.rows() != Eig_Modes.rows())
    {
        VolumesN.resize(Eig_Modes.rows(), 1);
        VolumesN.col(0).segment(0, Volumes.rows()) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() + 1, Volumes.rows() * 2) = Volumes;
        VolumesN.col(0).segment(Volumes.rows() * 2 + 1, Volumes.rows() * 3) = Volumes;
    }

    fr = Eig_Modes.transpose() * (f.cwiseProduct(VolumesN));
    return fr;
}
#endif
