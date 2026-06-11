#pragma once
#include <Eigen/Dense>

struct UKF_Weights
{
    Eigen::Matrix<double, 9, 1> Wm;
    Eigen::Matrix<double, 9, 1> Wc;
};

UKF_Weights UKF_compute_weights(int n, double alpha, double beta, double kappa);