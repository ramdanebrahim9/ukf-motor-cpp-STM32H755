#include "UKF_sigma_points.hpp"

Eigen::Matrix<double, 9, 4> UKF_compute_sigma_points(const Eigen::Vector4d& x,
                                                     const Eigen::Matrix4d& P,
                                                     double alpha,
                                                     double kappa)
{
    const int n = 4;
    const double lambda_ = alpha * alpha * (n + kappa) - n;

    Eigen::Matrix4d A = (lambda_ + n) * P;
    Eigen::LLT<Eigen::Matrix4d> llt(A);
    Eigen::Matrix4d U = llt.matrixL();

    Eigen::Matrix<double, 9, 4> sigmas;
    sigmas.row(0) = x.transpose();

    for (int k = 0; k < n; k++)
    {
        sigmas.row(k + 1) = (x + U.col(k)).transpose();
        sigmas.row(n + k + 1) = (x - U.col(k)).transpose();
    }

    return sigmas;
}