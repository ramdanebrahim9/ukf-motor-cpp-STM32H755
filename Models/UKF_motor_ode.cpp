#include "UKF_motor_ode.hpp"
#include <cmath>

Eigen::Vector2d UKF_motor_ode(const Eigen::Vector2d& x, double V)
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

    double omega = x(0);
    double i = x(1);

    // ── friction model ────────────────────────────────────────────────
    double sgn = tanh(omega / eps_s);
    double fric = sgn * (tcJ + (tsJ - tcJ) * exp(-pow(omega / ws, 2))) + BJ * omega;

    // ── ODEs ──────────────────────────────────────────────────────────
    double domega = (Ke / J) * i - fric;
    double di = (1.0 / L) * (V - R * i - Ke * omega);

    Eigen::Vector2d dx;
    dx(0) = domega;
    dx(1) = di;

    return dx;
}