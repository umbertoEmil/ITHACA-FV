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
\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the ITHACAregularization class, it contains the implementation of
/// several methods for regularization.

#include "ITHACAregularization.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace ITHACAregularization
{

Eigen::VectorXd  TSVD(Eigen::MatrixXd A,
                      Eigen::MatrixXd b, int filter)
{
    M_Assert(b.cols() == 1, "The b input in TSVD must have only one column");
    M_Assert(filter <= A.cols(), "Filter values too high");
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd x = Eigen::VectorXd::Zero(V.rows());

    for (label i = 0; i < filter; i++)
    {
        x += (U.col(i).dot(b.col(0)) / svd.singularValues()(i)) * (V.col(i));
    }

    return x;
}

Eigen::VectorXd  TSVD(Eigen::MatrixXd A,
                      Eigen::MatrixXd b, double noiseVariance, word parameterMethod)
{
    int filter;

    if (parameterMethod == "DP")
    {
        Info << "\nRegularization parameter selected by Discrepancy principle" << endl;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd U = svd.matrixU();
        Eigen::VectorXd bVect = b.col(0);
        double min;

        for (int col = 0; col < U.cols() - 1; col++)
        {
            double f = 0;

            for (int i = col + 1; i < U.cols(); i++)
            {
                Eigen::VectorXd tempU = U.col(i);
                f += tempU.col(i).dot(bVect) * tempU.col(i).dot(bVect);
            }

            f += 2 * noiseVariance * col;

            if (col == 0)
            {
                min = f;
                filter = col + 1;
            }
            else if (min > f)
            {
                min = f;
                filter = col + 1;
            }

            Info << "debug : f = " << f << endl;
            Info << "debug : min = " << min << endl;
            Info << "debug : k = " << filter << endl;
        }
    }
    else if (parameterMethod == "UPRE")
    {
        Info << "\nRegularization parameter selected by Discrepancy principle" << endl;
    }
    else
    {
        Info << "Regularization parameter selection methods available are:" << endl
             << "DP, UPRE" << endl;
        exit(1);
    }

    return ITHACAregularization::TSVD(A, b, filter);
}

Eigen::VectorXd  Tikhonov(Eigen::MatrixXd A,
                          Eigen::MatrixXd b, double regularizationParameter)
{
    M_Assert(b.cols() == 1, "The b input in Tikhonov regularization must have only one column");
    M_Assert(regularizationParameter >= 0, "The Tikhonov regularization parameter must be positive");
    Info << "RegularizationParameter = " << regularizationParameter << endl;

    Eigen::MatrixXd AtA = A.transpose() * A;
    Eigen::MatrixXd LHS = AtA + regularizationParameter * regularizationParameter * 
        Eigen::MatrixXd::Identity(AtA.rows(), AtA.cols());
    Eigen::MatrixXd RHS = A.transpose() * b;
    List<Eigen::MatrixXd> linSys;
    linSys.resize(2);
    linSys[0] = LHS;
    linSys[1] = RHS;
    Eigen::VectorXd x = linSys[0].fullPivLu().solve(linSys[1]);
    return x;
}

Eigen::VectorXd  Tikhonov(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd V,
                          Eigen::MatrixXd b, double regularizationParameter)
{
    M_Assert(b.cols() == 1, "The b input in Tikhonov regularization must have only one column");
    M_Assert(regularizationParameter >= 0, "The Tikhonov regularization parameter must be positive");
    Info << "RegularizationParameter = " << regularizationParameter << endl;

    label NsingVal = s.size();
    //Eigen::VectorXd beta = U.leftCols(NsingVal).transpose() * b;
    M_Assert(NsingVal == b.size(), 
            "b should have same size of the singular values vector");
    Eigen::VectorXd beta = U.leftCols(NsingVal).transpose() * b;
    Eigen::VectorXd zeta = s.cwiseProduct(beta);
    Eigen::VectorXd ss = s.cwiseProduct(s);
    Eigen::VectorXd temp = ss.array() + regularizationParameter * regularizationParameter;
    for(int i = 0; i < temp.size(); i++)
    {
        temp(i) = 1 / temp(i);
    }
    temp = zeta.cwiseProduct(temp); 
    Eigen::VectorXd x = V.leftCols(NsingVal) * temp;

    return x;
}

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd A, Eigen::MatrixXd b, label Nsteps, bool reorth, Eigen::VectorXd s)
{
    Info << "Using conjugate gradient method for the regularization of the linear system\n"
        << "Performing " << Nsteps 
        << " iterations applied implicitly to the normal equations A'*A*x = A'*b" 
        << endl << endl;
    M_Assert(b.cols() == 1, "The b input in conjugateGradient regularization must have only one column");
    M_Assert(Nsteps > 0, "conjugateGradient number of steps must be positive");
    if(reorth == 1)
    {
        M_Assert(s(0) > 0, "Singular values vector is needed for reorthogonalization");
        Info << "reorthogonalization not yet implemented, EXITING" << endl;
        exit(56);
    }

    //Initialization
    label Arows = A.rows();
    label Acols = A.cols();
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(Acols, Nsteps);
    Eigen::VectorXd eta = Eigen::VectorXd::Zero(Nsteps);
    Eigen::VectorXd rho = eta;

    //Setup for CG iteration
    Eigen::VectorXd x = Eigen::VectorXd::Zero(Acols);
    Eigen::VectorXd d = A.transpose() * b;
    Eigen::VectorXd r = b;
    scalar normr2 = d.squaredNorm();

    //Iterate
    for(label j = 0; j < Nsteps; j++)
    {
        //Update x and r vectors
        Eigen::VectorXd Ad = A * d;
        scalar alpha = normr2 / Ad.squaredNorm();
        x += alpha * d;
        r -= alpha * Ad;
        Eigen::VectorXd s = A.transpose() * r;

        //Update d vector
        scalar normr2_new = s.squaredNorm();
        scalar beta = normr2_new / normr2;
        normr2 = normr2_new;
        d = s + beta * d;
        X.col(j) = x;
        rho(j) = r.norm();
        eta(j) = x.norm();
    }

    word folder = "./ITHACAoutput/regularization/conjugateGradient";
    ITHACAstream::exportMatrix(X, "X", "eigen", folder);
    ITHACAstream::exportMatrix(rho, "residualNorm", "eigen", folder);
    ITHACAstream::exportMatrix(eta, "solutionNorm", "eigen", folder);

    return x;
}

Eigen::VectorXd PCGLS(Eigen::MatrixXd A, Eigen::MatrixXd L, Eigen::MatrixXd W, 
        Eigen::MatrixXd b, label Nsteps, bool reorth, Eigen::VectorXd singVal)
{
    scalar fudge_thr = 1e-4; // The fudge threshold is used to prevent filter factors from exploding
    // Initialization
    Info << "Using preconditioned conjugate gradient method\n" <<
        "for the regularization of the linear system\n"
        << "Performing " << Nsteps 
        << " iterations applied implicitly to the normal equations A'*A*x = A'*b" 
        << endl << endl;
    M_Assert(b.cols() == 1, "The b input in conjugateGradient regularization must have only one column");
    M_Assert(Nsteps > 0, "conjugateGradient number of steps must be positive");
    if(reorth == 1)
    {
        M_Assert(singVal(0) > 0, 
                "Singular values vector is needed for reorthogonalization");
        Info << "reorthogonalization not yet implemented, EXITING" << endl;
        exit(56);
    }
    label Arows = A.rows();
    label Acols = A.cols();
    label p = L.rows();
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(Acols, Nsteps);

    // Prepare for computations with L_p
    Eigen::MatrixXd S = pseudoInverse(A * W);
    Eigen::MatrixXd T = S * A;
    Eigen::VectorXd x = W * (S * b);
    std::cout << "debug : x = \n" << x << std::endl;

    // Prepare for CG iteration
    Eigen::VectorXd r = b - A * x;
    Info << "debug 1" << endl;
    Eigen::VectorXd s = A.transpose() * r;
    Info << "debug 2" << endl;
    Eigen::VectorXd q1 = ltsolve(L,s);
    Info << "q1.size() = " << q1.size() << endl;
    Eigen::VectorXd q = lsolve(L,q1, W, T);
    Info << "debug 4" << endl;
    Eigen::VectorXd z = q;
    Info << "debug 5" << endl;
    scalar dq = s.dot(q);
    Info << "debug 6" << endl;

    //Iterate
    for(label j = 0; j < Nsteps; j++)
    {
        Info << "Iteration " << j << endl;
        // Update x and r vectors; compute q1.
        Eigen::VectorXd Az  = A * z;
        scalar alpha = dq / Az.squaredNorm();
        x += alpha * z;
        r -= alpha * Az;
        s = A.transpose() * r;
        q1 = ltsolve(L,s);

        // Update z vector
        q = lsolve(L,q1, W, T);
        scalar dq2 = s.dot(q);
        scalar beta = dq2 / dq;
        dq = dq2;
        z = q + beta * z;
    }
    return x;
}

scalar GCV(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b, scalar regPar0,
        word method)
{
    M_Assert(method == "Tikhonov", "GCV is only implemented for Tikhonov regularization");
    label npoints = 1000;   //Number of points on the curve
    scalar smin_ratio = 16 * 2.2204e-16; //Smallest regularization parameter
    
    //Initialization
    label m = U.rows();
    label n = U.cols();
    Info << "debug : U = " << m << " x " << n << endl;
    label NsingVal = s.size();
    Eigen::VectorXd beta = U.transpose() * b;
    scalar beta2 = b.squaredNorm() - beta.squaredNorm();

    Eigen::VectorXd reg_param = Eigen::VectorXd::Zero(npoints); //Vector of regularization parameters
    Eigen::VectorXd G = reg_param;
    Eigen::VectorXd s2 = s.cwiseProduct(s);
    reg_param(npoints - 1) = std::max(s(NsingVal - 1), s(0) * smin_ratio);
    scalar ratio = std::pow(s(0) / reg_param(npoints - 1), 1.0 / (npoints * 1.0 - 1));
    for(int i = npoints - 2; i >= 0; i--)
    {
        reg_param(i) = ratio * reg_param(i+1);
    }
    scalar delta0 = 0; //Intrinsic residual
    if (m > n && beta2 > 0)
    {
        delta0 = beta2;
    }

    Eigen::VectorXd regPar(1);
    regPar(0) = regPar0;
    GCVfun_properties GCVprop;
    GCVprop.s2 = s2;
    GCVprop.beta = beta;
    GCVprop.delta0 = delta0;
    GCVprop.MnN = m - n;


    void* p = reinterpret_cast<void*>(&GCVprop);

    Eigen::VectorXd* grad_out;
    for(int i = 0; i < npoints; i++)
    {
        Eigen::VectorXd regI(1);
        regI(0) = reg_param(i);
        G(i) = GCVfun(regI, grad_out, p);
    }
    word folder = "./ITHACAoutput/GCV/";
    ITHACAstream::exportMatrix(G, "G", "eigen", folder);
    ITHACAstream::exportMatrix(reg_param, "regParam", "eigen", folder);


    optim::algo_settings_t settings;
    settings.vals_bound = true;
    Eigen::VectorXd regI(1);
    regI(0) = reg_param(npoints - 1);
    settings.lower_bounds = regI;
    regI(0) = reg_param(0);
    settings.upper_bounds = regI;

    std::cout << "debug : regPar = " << regPar << std::endl;

    bool success = optim::nm(regPar, GCVfun, p, settings);
    Info << "debug 1" << endl;
    if (success) {
        std::cout << "de: GCV minimum found.\n";
    } else {
        std::cout << "de: GCV minimum NOT found.\n Exiting" << std::endl;
        exit(13);
    }
    Info << "Regularization parameter = " << regPar(0) << endl;
    return regPar(0);
}

double
GCVfun(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out,
          void* opt_data)
{
    GCVfun_properties* objfn_data = reinterpret_cast<GCVfun_properties*>
                                    (opt_data);
    // Note: f = 1 - filter-factors.
    M_Assert(vals_inp.size() == 1, "Input vector to GCVfun must have size 1");
    scalar lambda = vals_inp(0);
    Eigen::VectorXd f(objfn_data->s2.size());
    for(int i = 0; i < f.size(); i++)
    {
        f(i) = lambda * lambda / (objfn_data->s2(i) + lambda * lambda);
    }
    Eigen::VectorXd temp = f.cwiseProduct(objfn_data->beta);
    return (temp.squaredNorm() + objfn_data->delta0) / 
        ((objfn_data->MnN + f.sum()) * (objfn_data->MnN + f.sum()));
}

void Lcurve(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b, word method, 
        label step)
{
    if(method != "Tikhonov")
    {
        Info << "Lcurve is only implemented for Tikhonov\n" <<
            "For CG method call Lcurve_CG" << endl << 
            "Exiting" << endl;
        exit(69);
    }

    label npoints;
    scalar smin_ratio;
    if(method == "Tikhonov")
    {
        npoints = 1000;   //Number of points on the curve
        smin_ratio = 16 * 2.2204e-16; //Smallest regularization parameter
    }
    
    //Initialization
    label m = U.rows();
    label n = U.cols();

    label NsingVal = s.size();
    Eigen::VectorXd beta = U.transpose() * b;
    beta = beta.head(NsingVal);
    scalar beta2 = b.squaredNorm() - beta.squaredNorm();
    Eigen::VectorXd xi = beta;
    for(int i = 0; i < xi.size(); i++)
    {
        xi(i) = xi(i) / s(i);
        if(!std::isfinite(xi(i)))
        {
            Info << "xi(" << i << ") in infinite" << endl;
            xi(i) = 0;
        }
    }

    Eigen::VectorXd eta = Eigen::VectorXd::Zero(npoints);
    Eigen::VectorXd rho = eta;
    Eigen::VectorXd reg_param = eta;
    Eigen::VectorXd s2 = s.cwiseProduct(s);
    reg_param(npoints - 1) = std::max(s(NsingVal - 1), s(0) * smin_ratio);
    Info << "debug : reg_param(npoints - 1) = " << reg_param(npoints - 1) << endl;
    Info << "debug : s(0) = " << s(0) << endl;
    double ratio = std::pow(s(0) / reg_param(npoints - 1), 1. / (npoints * 1.0 - 1.));
    Info << "debug : ratio = " << ratio << endl;
    for(int i = npoints - 2; i >= 0; i--)
    {
        reg_param(i) = ratio * reg_param(i + 1);
    }
    for(int i = 0; i < npoints; i++)
    {
        Eigen::VectorXd f(NsingVal);
        for(int j = 0; j < NsingVal; j++)
        {
            f(j) = s2(j) / (s2(j) + reg_param(i) * reg_param(i));
        }
        eta(i) = (f.cwiseProduct(xi)).norm();
        f = f - Eigen::VectorXd::Ones(f.size());
        rho(i) = (f.cwiseProduct(beta)).norm();
    }
    word folder = "./ITHACAoutput/Lcurve/";
    word etaout = "eta";
    etaout.append(std::to_string(step));
    word rhoout = "rho";
    rhoout.append(std::to_string(step));
    ITHACAstream::exportMatrix(eta, etaout, "eigen", folder);
    ITHACAstream::exportMatrix(rho, rhoout, "eigen", folder);
    ITHACAstream::exportMatrix(reg_param, "regParam", "eigen", folder);
    Info << "Lcurve done" << endl;
}

void Lcurve_CG(Eigen::MatrixXd A, Eigen::MatrixXd b)
{
    label npoints = 100;
    Eigen::VectorXd CGsteps(npoints);
    Eigen::VectorXd xNorm(npoints);
    Eigen::VectorXd AxbNorm(npoints);
    for(label k = 0; k < npoints; k++)
    {
        Eigen::VectorXd x = conjugateGradient(A, b, k + 1);
        xNorm(k) = x.norm();
        AxbNorm(k) = (A.transpose() * A * x - A.transpose() * b).norm();
    }
    word folder = "./ITHACAoutput/regularization/conjugateGradient/Lcurve";
    ITHACAstream::exportMatrix(xNorm, "xNorm", "eigen", folder);
    ITHACAstream::exportMatrix(AxbNorm, "AxbNorm", "eigen", folder);
    ITHACAstream::exportMatrix(CGsteps, "CGsteps", "eigen", folder);
}

void Picard(Eigen::MatrixXd U, Eigen::VectorXd s, Eigen::MatrixXd b)
{
    Info << "\nPreparing data for Picard plot" << endl;
    int NsingVal = s.size();
    Eigen::VectorXd beta = U.leftCols(NsingVal).transpose() * b;
    beta = beta.cwiseAbs();
    Eigen::VectorXd eta = Eigen::VectorXd::Zero(NsingVal);
    for(int i = 0; i < NsingVal; i++)
    {
        eta(i) = beta(i)/s(i);
    }
    word folder = "./ITHACAoutput/Picard/";
    ITHACAstream::exportMatrix(eta, "eta", "eigen", folder);
    ITHACAstream::exportMatrix(s, "singVal", "eigen", folder);
    ITHACAstream::exportMatrix(beta, "beta", "eigen", folder);
    Info << "Picard data ready\n" << endl;
}

Eigen::MatrixXd get_L(label n, label d, Eigen::Ref<Eigen::MatrixXd> W)
{
    M_Assert(d>=0, "Order d must be nonnegative");

    // Zero'th derivative.
    if (d==0)
    {
        Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n,n);
        W = Eigen::MatrixXd::Zero(n,1);
        return L;
    }
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n-1,n);
    if (d==1)
    {
        for(label i = 0; i < n -1; i++)
        {
            for(label j = 0; j < n; j++)
            {
                if(i == j)
                {
                    L(i,j) = -1;
                }
                else if(j == i + 1)
                {
                    L(i,j) = 1;
                }
            }
        }
    }
    else
    {
        Info << "ITHACAregularization::get_L not yet implemented for d > 1" << endl << 
            "Exiting" << endl;
        exit(78);
    }
    
    M_Assert(W.rows() == n && W.cols() == d, "W has wrong size");
    W = Eigen::MatrixXd::Zero(n,d);
    W.col(0) = Eigen::MatrixXd::Ones(n,1);
    return L;
}

Eigen::MatrixXd pseudoInverse(Eigen::MatrixXd a, double epsilon)
{
    // For a non-square matrix
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = 
        epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    return svd.matrixV() * 
        (svd.singularValues().array().abs() > tolerance).select(
                svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * 
        svd.matrixU().adjoint();
}

Eigen::VectorXd ltsolve(Eigen::MatrixXd L, Eigen::VectorXd y, Eigen::MatrixXd W, Eigen::MatrixXd T)
{
    Eigen::VectorXd x;
    if(L.cols() == L.rows()) // Square L
    {
        x = L.transpose().fullPivLu().solve(y);
        return x;
    }

    if(W.size() > 0)
    {
        y = y.head(L.rows()) - T.leftCols(L.rows()).transpose() * (W.transpose() * y);
        std::cout << "debug ltsolve : y = \n" << y << std::endl;
        Info << "y.size() = " << y.size() << endl;
        
    }
    
    x = (L.leftCols(L.rows()).transpose()).fullPivLu().solve(y.head(L.rows()));
    return x;
}

Eigen::VectorXd lsolve(Eigen::MatrixXd L, Eigen::VectorXd y, Eigen::MatrixXd W, Eigen::MatrixXd T)
{
    Eigen::VectorXd x;
    if(L.cols() == L.rows()) // Square L
    {
        x = L.transpose().fullPivLu().solve(y);
        return x;
    }
    Info << "y.size() = " << y.size() << endl;
    Info << "L = " << L.rows() << " x " << L.cols() << endl;

    x = (L.leftCols(L.rows()).transpose()).fullPivLu().solve(y);
    Info << "debug lsolve 1" << endl;
    Info << "T = " << T.rows() << " x " << T.cols() << endl;
    std::cout << "W * (T.leftCols(L.rows()) * x) = \n" << W * (T.leftCols(L.rows()) * x) << std::endl;
    Info << "x.size() = " << x.size() << endl;
    Info << "W * (T.leftCols(L.rows()) * x) = " << (W * (T.leftCols(L.rows()) * x)).size() << endl;
    Eigen::VectorXd xTemp(x.size() + (L.cols() - L.rows()));
    xTemp << x, Eigen::VectorXd::Zero(L.cols() - L.rows());
    x = xTemp - W * (T.leftCols(L.rows()) * x); 
    return x;
}

List<Eigen::MatrixXd> precondition(List<Eigen::MatrixXd> linSys, word method, int nSVD)
{
    M_Assert(linSys.size() == 2, "The linear system has wrong size");
    M_Assert(method != "None", "Specify the preconditioner to use");

    List<Eigen::MatrixXd> precLS = linSys;
    Eigen::MatrixXd precA;
    if(method == "Jacobi")
    {
        Eigen::VectorXd temp = linSys[0].diagonal();
        for(int i = 0; i < temp.size(); i++)
        {
            temp(i) = 1.0 / temp(i);
        }
        precA = temp.asDiagonal();
    }
    else if(method == "SVD")
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::VectorXd singVal = svd.singularValues();
        M_Assert(nSVD > 0, "Wrong input value to SVD preconditioning method");
        for(int i = 0; i < singVal.size(); i++)
        {
            if(i >= nSVD)
            {
                singVal(i) = 1;
            }
        }
        precA = svd.matrixU() * singVal.asDiagonal() * svd.matrixV().transpose();
        precA = precA.inverse();
    }
    precLS[0] = precA * linSys[0];
    precLS[1] = precA * linSys[1];

    word outputFolder = "./ITHACAoutput/regularization/preconditioning";
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(linSys[0],
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singVal = svd.singularValues();
    ITHACAstream::exportMatrix(singVal, "singularValues", "eigen", outputFolder);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_precond(precLS[0],
                                          Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd singVal_precond = svd_precond.singularValues();
    ITHACAstream::exportMatrix(singVal_precond, "singularValues_precond", "eigen", outputFolder);
    ITHACAstream::exportMatrix(precLS[0], "A_prec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(precLS[1], "b_prec", "eigen", outputFolder);
    ITHACAstream::exportMatrix(linSys[0], "A", "eigen", outputFolder);
    ITHACAstream::exportMatrix(linSys[1], "b", "eigen", outputFolder);
    return precLS;

        
}


}
