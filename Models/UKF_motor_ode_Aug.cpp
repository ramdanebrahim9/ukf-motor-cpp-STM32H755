#include "UKF_motor_ode_Aug.hpp"
#include <cmath>

Eigen::Vector4d UKF_motor_ode_Aug(const Eigen::Vector4d& x, double V,
                                  double i_rk4_k, double omega_rk4_k)
{
    // ── motor parameters ──────────────────────────────────────────────
    static constexpr double R = 1.438907;
    static constexpr double L = 0.415e-3;
    static constexpr double Ke = 0.137961;
    static constexpr double J = 1.263076e-4;
    static constexpr double BJ = 6.895330;
    static constexpr double tcJ = 117.827967;
    static constexpr double tsJ = 104.736000;
    static constexpr double ws = 5.0;
    static constexpr double eps_s = 0.5;
    static constexpr double tau_d = 0.1;

    double theta = x(0); // estimated
    double omega = x(1); // estimated
    double i = x(2);     // estimated
    double d = x(3);     // estimated

    // ── mechanical side — i_rk4 injected ─────────────────────────────
    double sgn = tanh(omega / eps_s);
    double tmp = omega / ws;
    double fric = sgn * (tcJ + (tsJ - tcJ) * exp(-(tmp * tmp))) + BJ * omega;
    double dtheta = omega;
    double domega = (Ke / J) * i_rk4_k - fric;

    // ── electrical side — omega_rk4 injected ─────────────────────────
    double di = (1.0 / L) * (V - R * i - Ke * omega_rk4_k + d);

    // ── disturbance — untouched ───────────────────────────────────────
    double dd = -d / tau_d;

    Eigen::Vector4d dx;
    dx(0) = dtheta;
    dx(1) = domega;
    dx(2) = di;
    dx(3) = dd;

    return dx;
}