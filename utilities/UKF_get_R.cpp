#include "UKF_get_R.hpp"

Eigen::Matrix3d UKF_get_R(double V, double omega_rk4, bool tick_fired)
{
    // ── current noise — adaptive on omega ─────────────────────────────
    static const double w_table[10] = {6.5, 10.71, 17.0, 25.0, 31.0, 40.69, 47.94, 54.4, 63.5, 70.10};
    static const double R_table[10] = {0.095 * 2, 0.11 * 2, 0.145 * 2, 0.16 * 2, 0.16 * 2,
                                       0.1675 * 2, 0.185 * 2, 0.16 * 2, 0.145 * 2, 0.13 * 2};

    double w_clamped = fmax(w_table[0], fmin(w_table[9], omega_rk4));
    double r_i = UKF_interp1(w_table, R_table, 10, w_clamped);

    // ── theta + omega noise — adaptive on V, tick-scheduled ───────────
    static const double V_bp[7] = {0.0, 6.5, 7.5, 8.5, 10.0, 11.0, 12.0};
    static const double k1[7] = {0.15, 0.15, 0.90, 0.90, 0.90, 1.50, 2.00}; // trusted omega std
    static const double k2[7] = {1.0, 1.0, 4.0, 4.0, 4.0, 5.0, 5.0};        // coasting omega std

    double V_c = fmax(V_bp[0], fmin(V_bp[6], fabs(V)));
    double k1_v = UKF_interp1(V_bp, k1, 7, V_c);
    double k2_v = UKF_interp1(V_bp, k2, 7, V_c);

    static constexpr double tick_rad = (2.0 * M_PI) / (64.0 * 13.7335);

    double r_theta, r_omega;
    if (tick_fired)
    {
        r_theta = 1e-6; // sqrt(1e-12)
        r_omega = k1_v; // trusted omega std
    }
    else
    {
        r_theta = tick_rad;
        r_omega = k2_v; // coasting omega std
    }

    // ── assemble 3x3 R — [theta, omega, i] ───────────────────────────
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
    R(0, 0) = r_theta * r_theta;
    R(1, 1) = r_omega * r_omega;
    R(2, 2) = r_i * r_i;

    return R;
}