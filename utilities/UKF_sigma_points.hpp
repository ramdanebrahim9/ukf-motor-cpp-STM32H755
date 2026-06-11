#pragma once
#include <Eigen/Dense>

Eigen::Matrix<double, 9, 4> UKF_compute_sigma_points(const Eigen::Vector4d& x,
                                                     const Eigen::Matrix4d& P,
                                                     double alpha,
                                                     double kappa);