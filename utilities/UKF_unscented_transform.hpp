#pragma once
#include <Eigen/Dense>

template <int N>
struct UKF_UTResult
{
    Eigen::Matrix<double, N, 1> x; // weighted mean
    Eigen::Matrix<double, N, N> P; // weighted covariance + noise
};

template <int N>
UKF_UTResult<N> UKF_unscented_transform(
    const Eigen::Matrix<double, 9, N>& sigmas,
    const Eigen::Matrix<double, 9, 1>& Wm,
    const Eigen::Matrix<double, 9, 1>& Wc,
    const Eigen::Matrix<double, N, N>& noise_cov)
{
    UKF_UTResult<N> out;

    // weighted mean
    out.x = sigmas.transpose() * Wm;

    // residuals — each row = sigma_i - mean
    Eigen::Matrix<double, 9, N> y = sigmas.rowwise() - out.x.transpose();

    // weighted covariance
    out.P = y.transpose() * Wc.asDiagonal() * y;

    // add noise
    out.P += noise_cov;

    return out;
}