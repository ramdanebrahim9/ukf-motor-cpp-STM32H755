#pragma once
#include <Eigen/Dense>

Eigen::Vector2d UKF_motor_ode(const Eigen::Vector2d& x, double V);