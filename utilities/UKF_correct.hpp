#pragma once
#include <Eigen/Dense>
#include "UKF_weights.hpp"

struct UKF_CorrectResult
{
    Eigen::Vector4d x_est;         // (4x1)
    Eigen::Matrix4d P_est;         // (4x4)
    Eigen::Matrix<double, 4, 3> K; // (4x3)
    Eigen::Vector3d innov;         // (3x1) — y in MATLAB
};

UKF_CorrectResult UKF_correct(const Eigen::Vector4d& x_pred,
                              const Eigen::Matrix4d& P_pred,
                              const Eigen::Matrix<double, 9, 4>& sigmas_f,
                              const Eigen::Vector3d& z,
                              const Eigen::Matrix3d& R,
                              const UKF_Weights& W);