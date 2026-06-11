#include "UKF_predict.hpp"
#include "UKF_sigma_points.hpp"
#include "UKF_rk4_step_Aug.hpp"
#include "UKF_unscented_transform.hpp"

UKF_PredictResult UKF_predict(const Eigen::Vector4d& x,
                              const Eigen::Matrix4d& P,
                              const Eigen::Matrix4d& Q,
                              const UKF_Weights& W,
                              double alpha,
                              double kappa,
                              double dt,
                              double V,
                              double i_rk4_k,
                              double omega_rk4_k)
{
    // ── generate sigma points ─────────────────────────────────────────
    Eigen::Matrix<double, 9, 4> sigmas = UKF_compute_sigma_points(x, P, alpha, kappa);

    // ── propagate each sigma point through f (rk4_step_Aug) ──────────
    Eigen::Matrix<double, 9, 4> sigmas_f;
    for (int i = 0; i < 9; i++)
    {
        Eigen::Vector4d sigma_i = sigmas.row(i).transpose();
        sigmas_f.row(i) = UKF_rk4_step_Aug(sigma_i, V, i_rk4_k, omega_rk4_k, dt).transpose();
    }

    // ── unscented transform → predicted mean and covariance ───────────
    UKF_UTResult<4> ut = UKF_unscented_transform<4>(sigmas_f, W.Wm, W.Wc, Q);

    UKF_PredictResult res;
    res.x_pred = ut.x;
    res.P_pred = ut.P;
    res.sigmas_f = sigmas_f;

    return res;
}