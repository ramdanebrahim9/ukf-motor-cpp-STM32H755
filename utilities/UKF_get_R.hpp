#pragma once
#include <Eigen/Dense>
#include <cmath>
#include "UKF_utils.hpp"

Eigen::Matrix3d UKF_get_R(double V, double omega_rk4, bool tick_fired);