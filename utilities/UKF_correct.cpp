#include "UKF_correct.hpp"
#include "UKF_unscented_transform.hpp"

UKF_CorrectResult UKF_correct(const Eigen::Vector4d& x_pred,
                              const Eigen::Matrix4d& P_pred,
                              const Eigen::Matrix<double, 9, 4>& sigmas_f,
                              const Eigen::Vector3d& z,
                              const Eigen::Matrix3d& R,
                              const UKF_Weights& W)
{
    // ── propagate sigmas through h — h(x) = [theta, omega, i] = x(0:2) ──
    Eigen::Matrix<double, 9, 3> sigmas_h;
    for (int i = 0; i < 9; i++)
        sigmas_h.row(i) = sigmas_f.row(i).head<3>();

    // ── unscented transform → predicted measurement mean and covariance ──
    UKF_UTResult<3> ut = UKF_unscented_transform<3>(sigmas_h, W.Wm, W.Wc, R);
    Eigen::Vector3d zp = ut.x; // predicted measurement
    Eigen::Matrix3d S = ut.P;  // innovation covariance

    // ── cross covariance Pxz ─────────────────────────────────────────────
    Eigen::Matrix<double, 4, 3> Pxz = Eigen::Matrix<double, 4, 3>::Zero();
    for (int i = 0; i < 9; i++)
    {
        Eigen::Vector4d dx = sigmas_f.row(i).transpose() - x_pred;
        Eigen::Vector3d dz = sigmas_h.row(i).transpose() - zp;
        Pxz += W.Wc(i) * (dx * dz.transpose());
    }

    // ── Kalman gain ───────────────────────────────────────────────────────
    Eigen::Matrix<double, 4, 3> K = Pxz * S.inverse();

    // ── innovation ────────────────────────────────────────────────────────
    Eigen::Vector3d y = z - zp;

    // ── update state and covariance ───────────────────────────────────────
    UKF_CorrectResult res;
    res.x_est = x_pred + K * y;
    res.P_est = P_pred - K * S * K.transpose();
    res.K = K;
    res.innov = y;

    return res;
}