#include "UKF_rk4_step_Aug.hpp"
#include "UKF_motor_ode_Aug.hpp"

Eigen::Vector4d UKF_rk4_step_Aug(const Eigen::Vector4d& x, double V,
                                 double i_rk4_k, double omega_rk4_k, double dt)
{
    Eigen::Vector4d k1 = UKF_motor_ode_Aug(x, V, i_rk4_k, omega_rk4_k);
    Eigen::Vector4d k2 = UKF_motor_ode_Aug(x + (dt / 2.0) * k1, V, i_rk4_k, omega_rk4_k);
    Eigen::Vector4d k3 = UKF_motor_ode_Aug(x + (dt / 2.0) * k2, V, i_rk4_k, omega_rk4_k);
    Eigen::Vector4d k4 = UKF_motor_ode_Aug(x + dt * k3, V, i_rk4_k, omega_rk4_k);

    return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}