#pragma once
#include <Eigen/Dense>
#include "UKF_weights.hpp"
#include "UKF_predict.hpp"
#include "UKF_correct.hpp"
#include "UKF_get_Q.hpp"
#include "UKF_get_R.hpp"
#include "UKF_rk4_step.hpp"

// ── constants ─────────────────────────────────────────────────────────────
static constexpr double UKF_TICKS_TO_RAD = (2.0 * M_PI) / (64.0 * 13.7335);
static constexpr double UKF_TS = 50e-6;
static constexpr double UKF_ALPHA = 1e-3;
static constexpr double UKF_BETA = 2.0;
static constexpr double UKF_KAPPA = 0.0;

// ── IIR low-pass: fc = 200 Hz, Ts = 50us ─────────────────────────────────
// alpha = exp(-2*pi*200*50e-6) ≈ 0.9391
static constexpr double UKF_IIR_ALPHA = 0.9391;

struct UKF_State
{
    Eigen::Vector4d x_est; // [theta, omega, i, d]
    Eigen::Matrix4d P_est;
};

// ── initial conditions (optional) ────────────────────────────────────────
struct UKF_InitParams
{
    double theta_0 = 0.0;
    double omega_0 = 0.0;
    double i_0 = 0.0;
    double d_0 = 0.0;
};

class UKF
{
  public:
    // ── init ──────────────────────────────────────────────────────────────
    static void init();                         // all initial conditions  zeros
    static void init(const UKF_InitParams& ic); // warm start

    // ── main step — call once per control loop ────────────────────────────
    // inputs : raw sample fields
    // returns: current best estimate
    static const UKF_State& step(double u,
                                 double i_meas,
                                 int32_t ticks,
                                 double omega_dirty);

    // ── getters ───────────────────────────────────────────────────────────
    static const UKF_State& getState()
    {
        return s_state;
    }
    static const Eigen::Vector2d& getRK4State()
    {
        return s_x_rk4;
    }

    // ── UKF state ─────────────────────────────────────────────────────────
    static UKF_State s_state;

    // ── open-loop RK4 model state [omega, i] ─────────────────────────────
    static Eigen::Vector2d s_x_rk4;

    // ── tick tracking ─────────────────────────────────────────────────────
    static int32_t s_last_ticks;

    // ── IIR filter state for mean_I (disturbance clamp) ──────────────────
    static double s_iir_mean_i;

    // ── weights — computed once at init ───────────────────────────────────
    static UKF_Weights s_W;

    // ── disturbance offset tables (from MATLAB characterization) ──────────
    static const double s_V_offset_table[12];
    static const double s_offset_table[12];

    // ── debug / test snapshots (written each step) ────────────────────────
    static bool s_tick_fired;
    static Eigen::Matrix4d s_Q_k;
    static Eigen::Matrix3d s_R_k;
    static Eigen::Vector3d s_z_k;
    static Eigen::Matrix<double, 9, 4> s_sigmas;
    static Eigen::Vector4d s_x_pred;
    static Eigen::Matrix4d s_P_pred;
    static Eigen::Vector3d s_innov;
    static Eigen::Matrix<double, 4, 3> s_K;

    // ── internal helpers ──────────────────────────────────────────────────
    static double iir_update(double x_new);
    static bool detect_tick(int32_t ticks);
    static double ticks_to_rad(int32_t ticks);
    static double compute_disturbance_clip(double V, double mean_i, double i_rk4);
};

template <typename Derived>
void printMatrix(const Eigen::MatrixBase<Derived>& mat, const char* name)
{
    printf("%s:\n\r", name);
    for (int i = 0; i < mat.rows(); i++)
    {
        for (int j = 0; j < mat.cols(); j++)
        {
            printf("%.10e ", (double)mat(i, j));
        }
        printf("\n\r");
    }
    printf("\n\r");
};