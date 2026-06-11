#pragma once
#include <Eigen/Dense>

Eigen::Vector4d UKF_motor_ode_Aug(const Eigen::Vector4d& x, double V,
                                  double i_rk4_k, double omega_rk4_k);