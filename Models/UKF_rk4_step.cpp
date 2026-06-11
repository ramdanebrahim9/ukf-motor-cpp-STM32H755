#include "UKF_rk4_step.hpp"
#include "UKF_motor_ode.hpp"

Eigen::Vector2d UKF_rk4_step(const Eigen::Vector2d& x, double V, double dt)
{
    Eigen::Vector2d k1 = UKF_motor_ode(x, V);
    Eigen::Vector2d k2 = UKF_motor_ode(x + (dt / 2.0) * k1, V);
    Eigen::Vector2d k3 = UKF_motor_ode(x + (dt / 2.0) * k2, V);
    Eigen::Vector2d k4 = UKF_motor_ode(x + dt * k3, V);

    return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}