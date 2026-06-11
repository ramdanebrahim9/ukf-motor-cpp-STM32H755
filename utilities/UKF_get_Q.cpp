#include "UKF_get_Q.hpp"

Eigen::Matrix4d UKF_get_Q(double V)
{
    static const double V_table[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.5, 7.5, 8.5, 9.0, 10.0};
    static const double Q_table[10] = {
        0.025 * 0.45, 0.025 * 0.45, 0.025 * 0.45, 0.025 * 0.45,
        0.022 * 0.45, 0.033 * 0.45, 0.030 * 0.45, 0.033 * 0.45,
        0.0255 * 0.45, 0.018 * 0.45};

    // clamp
    double V_clamped = fmax(V_table[0], fmin(V_table[9], V));

    // interpolate
    double q_i = UKF_interp1(V_table, Q_table, 10, V_clamped);

    double q_d = 0.01;

    Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
    Q(0, 0) = 1e-8;      // theta — variance directly
    Q(1, 1) = 1e-4;      // omega — variance directly
    Q(2, 2) = q_i * q_i; // i     — adaptive, squared std
    Q(3, 3) = q_d * q_d; // d     — squared std

    return Q;
}