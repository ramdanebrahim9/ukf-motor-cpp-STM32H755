#include "UKF.hpp"
#include <cmath>

// ── static member definitions ─────────────────────────────────────────────
UKF_State UKF::s_state;
Eigen::Vector2d UKF::s_x_rk4;
int32_t UKF::s_last_ticks = 0;
double UKF::s_iir_mean_i = 0.0;
UKF_Weights UKF::s_W;

bool UKF::s_tick_fired = false;
Eigen::Matrix4d UKF::s_Q_k;
Eigen::Matrix3d UKF::s_R_k;
Eigen::Vector3d UKF::s_z_k;
Eigen::Matrix<double, 9, 4> UKF::s_sigmas;
Eigen::Vector4d UKF::s_x_pred;
Eigen::Matrix4d UKF::s_P_pred;
Eigen::Vector3d UKF::s_innov;
Eigen::Matrix<double, 4, 3> UKF::s_K;

const double UKF::s_V_offset_table[12] = {
    1.193, 1.852, 2.881, 3.946, 5.077, 6.286,
    7.459, 8.404, 9.207, 10.034, 10.927, 11.902};
const double UKF::s_offset_table[12] = {
    0.0127, -0.0024, -0.0181, -0.0271, -0.0269, -0.0331,
    -0.0466, -0.0257, -0.0277, -0.0743, -0.0908, -0.1680};

void UKF::init()
{
    init(UKF_InitParams{}); // delegate — all defaults are 0.0
}

void UKF::init(const UKF_InitParams& ic)
{
    s_W = UKF_compute_weights(4, UKF_ALPHA, UKF_BETA, UKF_KAPPA);

    s_state.x_est(0) = ic.theta_0;
    s_state.x_est(1) = ic.omega_0;
    s_state.x_est(2) = ic.i_0;
    s_state.x_est(3) = ic.d_0;

    s_state.P_est = Eigen::Matrix4d::Zero();
    s_state.P_est(0, 0) = 0.01 * 0.01;
    s_state.P_est(1, 1) = 0.5 * 0.5;
    s_state.P_est(2, 2) = 0.1 * 0.1;
    s_state.P_est(3, 3) = 0.08 * 0.08;

    s_x_rk4(0) = ic.omega_0; // RK4 shares same initial omega and i
    s_x_rk4(1) = ic.i_0;

    s_last_ticks = 0;
    s_iir_mean_i = ic.i_0; // seed IIR from measured current, not 0
}

// ── iir_update ────────────────────────────────────────────────────────────
double UKF::iir_update(double x_new)
{
    s_iir_mean_i = UKF_IIR_ALPHA * s_iir_mean_i + (1.0 - UKF_IIR_ALPHA) * x_new;
    return s_iir_mean_i;
}

// ── detect_tick ───────────────────────────────────────────────────────────
bool UKF::detect_tick(int32_t ticks)
{
    bool fired = (ticks != s_last_ticks);
    s_last_ticks = ticks;
    return fired;
}

// ── ticks_to_rad ──────────────────────────────────────────────────────────
double UKF::ticks_to_rad(int32_t ticks)
{
    return static_cast<double>(ticks) * UKF_TICKS_TO_RAD;
}

// ── compute_disturbance_clip ──────────────────────────────────────────────
double UKF::compute_disturbance_clip(double V, double mean_i, double i_rk4)
{
    double V_c = fmax(s_V_offset_table[0], fmin(s_V_offset_table[11], fabs(V)));
    double offset_k = fabs(UKF_interp1(s_V_offset_table, s_offset_table, 12, V_c));

    double thresh_k;
    if (offset_k >= 0.06)
        thresh_k = 0.035;
    else if (offset_k >= 0.039)
        thresh_k = 0.042;
    else
        thresh_k = 0.0615;

    double drift = mean_i - i_rk4;
    return (fabs(drift) > thresh_k) ? fabs(drift) : 0.0;
}

// ── step ─────────────────────── ──────────────────────────────────────────
const UKF_State& UKF::step(double u,
                           double i_meas,
                           int32_t ticks,
                           double omega_dirty)
{
    // ── 1. tick detection and theta ───────────────────────────────────────
    s_tick_fired = detect_tick(ticks);
    double theta_k = ticks_to_rad(ticks);

    // ── 2. IIR filter on current ──────────────────────────────────────────
    double mean_i_k = iir_update(i_meas);

    // ── 3. open-loop RK4 step — produces i_rk4_k, omega_rk4_k ───────────
    double omega_rk4_k = s_x_rk4(0);
    double i_rk4_k = s_x_rk4(1);
    s_x_rk4 = UKF_rk4_step(s_x_rk4, u, UKF_TS); // advance for next step

    // ── 4. adaptive Q and R ───────────────────────────────────────────────
    s_Q_k = UKF_get_Q(u);
    s_R_k = UKF_get_R(u, omega_rk4_k, s_tick_fired);

    // ── 5. measurement vector z = [theta, omega_dirty, i_meas] ───────────
    // omega_dirty: use rk4 on tick, hold UKF estimate between ticks
    double omega_z = s_tick_fired ? omega_rk4_k : s_state.x_est(1);
    s_z_k(0) = theta_k;
    s_z_k(1) = omega_z;
    s_z_k(2) = i_meas;

    // ── 6. UKF predict ────────────────────────────────────────────────────
    UKF_PredictResult pred = UKF_predict(
        s_state.x_est, s_state.P_est,
        s_Q_k, s_W,
        UKF_ALPHA, UKF_KAPPA,
        UKF_TS, u,
        i_rk4_k, omega_rk4_k);

    s_x_pred = pred.x_pred;
    s_P_pred = pred.P_pred;
    s_sigmas = pred.sigmas_f;

    // ── 7. UKF correct ────────────────────────────────────────────────────
    UKF_CorrectResult corr = UKF_correct(
        pred.x_pred, pred.P_pred,
        pred.sigmas_f,
        s_z_k, s_R_k, s_W);

    s_innov = corr.innov;
    s_K = corr.K;

    s_state.x_est = corr.x_est;
    s_state.P_est = corr.P_est;

    // ── 8. disturbance clamp ──────────────────────────────────────────────
    double d_clip = compute_disturbance_clip(u, mean_i_k, i_rk4_k);
    double d = s_state.x_est(3);
    s_state.x_est(3) = fmax(fmin(d, d_clip), -d_clip);

    return s_state;
}