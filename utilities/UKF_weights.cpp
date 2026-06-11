#include "UKF_weights.hpp"

UKF_Weights UKF_compute_weights(int n, double alpha, double beta, double kappa)
{
    const double lambda_ = alpha * alpha * (n + kappa) - n;
    const double c = 0.5 / (n + lambda_);

    UKF_Weights W;
    W.Wm = Eigen::Matrix<double, 9, 1>::Constant(c);
    W.Wc = Eigen::Matrix<double, 9, 1>::Constant(c);

    W.Wm(0) = lambda_ / (n + lambda_);
    W.Wc(0) = lambda_ / (n + lambda_) + (1.0 - alpha * alpha + beta);

    return W;
}