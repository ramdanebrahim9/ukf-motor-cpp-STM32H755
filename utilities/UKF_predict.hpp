#pragma once
#include <Eigen/Dense>
#include "UKF_weights.hpp"

struct UKF_PredictResult
{
    Eigen::Vector4d x_pred;               // (4x1)
    Eigen::Matrix4d P_pred;               // (4x4)
    Eigen::Matrix<double, 9, 4> sigmas_f; // (9x4)
};

UKF_PredictResult UKF_predict(const Eigen::Vector4d& x,
                              const Eigen::Matrix4d& P,
                              const Eigen::Matrix4d& Q,
                              const UKF_Weights& W,
                              double alpha,
                              double kappa,
                              double dt,
                              double V,
                              double i_rk4_k,
                              double omega_rk4_k);