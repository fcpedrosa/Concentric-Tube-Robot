#pragma once

// Internal math utilities of the CTR kinematics library.
// This header lives under detail/ — it is installed (public headers include
// it) but its contents are NOT part of the stable public API.

#include <cmath>
#include <numbers>
#include <blaze/Math.h>

namespace ctr::math
{

/**
 * @brief Converts an angle in degrees to radians.
 *
 * @param degree Angle in degrees.
 * @return The corresponding angle in radians.
 */
inline constexpr double deg2Rad(double degree)
{
    constexpr double pi_180 = std::numbers::pi / 180.00;
    return degree * pi_180;
}

/**
 * @brief Wraps an angle to the interval (-π, π] by subtracting an exact
 *        multiple of 2π.
 *
 * Because the reduction is an exact 2π-congruence (up to floating-point
 * rounding), any 2π-periodic function of the angle — in particular the CTR
 * forward kinematics as a function of the tube rotations α — is invariant
 * under this wrap.
 *
 * @param angle Angle in radians.
 * @return The congruent angle in (-π, π].
 */
inline double wrapToPi(double angle) noexcept
{
    return std::remainder(angle, 2.0 * std::numbers::pi);
}

/**
 * @brief Converts an Euler-angle representation of 3D rotations into the corresponding quaternion.
 * The quaternion is passed by reference and is modified by the function.
 *
 * @param heading  The heading angle [rad].
 * @param attitude The attitude angle [rad].
 * @param bank     The bank angle [rad].
 * @param h        Output: the corresponding quaternion [w, x, y, z].
 */
inline void euler2Quaternion(const double heading, const double attitude, const double bank,
                             blaze::StaticVector<double, 4UL> &h)
{
    const double theta(0.50 * heading), phi(0.50 * attitude), psi(0.50 * bank);
    const double c1(std::cos(theta)), s1(std::sin(theta)), c2(std::cos(phi)), s2(std::sin(phi)), c3(std::cos(psi)),
        s3(std::sin(psi));

    const double c1c2 = c1 * c2;
    const double s1s2 = s1 * s2;

    h[0UL] = c1c2 * c3 - s1s2 * s3;
    h[1UL] = c1c2 * s3 + s1s2 * c3;
    h[2UL] = s1 * c2 * c3 + c1 * s2 * s3;
    h[3UL] = c1 * s2 * c3 - s1 * c2 * s3;
}

/**
 * @brief Converts a 3D rotation expressed as a quaternion into the corresponding rotation matrix in SO(3).
 *
 * Accepts any dense-vector expression (e.g. a blaze::subvector view into the
 * ODE state) so callers do not materialize a temporary StaticVector.
 * The scale factor 2/||h||² makes the extracted matrix self-normalizing, so
 * modest quaternion drift during integration is absorbed here.
 *
 * @param h The input 4-dimensional quaternion [w, x, y, z] (any dense vector expression).
 * @param R Output: the rotation matrix in SO(3).
 */
template <typename VT>
inline void getSO3(const blaze::DenseVector<VT, blaze::columnVector> &h_,
                   blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> &R)
{
    const auto &h = *h_;
    const double scale = 2.00 / blaze::sqrNorm(h);

    R(0UL, 0UL) = 1.00 + scale * (-h[2UL] * h[2UL] - h[3UL] * h[3UL]);
    R(0UL, 1UL) = scale * (h[1UL] * h[2UL] - h[3UL] * h[0UL]);
    R(0UL, 2UL) = scale * (h[1UL] * h[3UL] + h[2UL] * h[0UL]);

    R(1UL, 0UL) = scale * (h[1UL] * h[2UL] + h[3UL] * h[0UL]);
    R(1UL, 1UL) = 1.00 + scale * (-h[1UL] * h[1UL] - h[3UL] * h[3UL]);
    R(1UL, 2UL) = scale * (h[2UL] * h[3UL] - h[1UL] * h[0UL]);

    R(2UL, 0UL) = scale * (h[1UL] * h[3UL] - h[2UL] * h[0UL]);
    R(2UL, 1UL) = scale * (h[2UL] * h[3UL] + h[1UL] * h[0UL]);
    R(2UL, 2UL) = 1.00 + scale * (-h[1UL] * h[1UL] - h[2UL] * h[2UL]);
}

/**
 * @brief Computes and returns the rotation matrix in SO(3) about the z-axis for an angle 'theta' in radians.
 *
 * @param theta The input angle in radians.
 * @return The corresponding rotation matrix about the z-axis.
 */
inline blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> rotz(const double theta)
{
    const double c(std::cos(theta)), s(std::sin(theta));
    return blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>{
        {c, -s, 0.00}, {s, c, 0.00}, {0.00, 0.00, 1.00}};
}

/**
 * @brief Computes the hat (skew-symmetric) operator applied to a 3-dimensional vector.
 *
 * @param v The input 3-dimensional vector.
 * @return The corresponding skew-symmetric matrix in so(3).
 */
inline blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> hatOperator(const blaze::StaticVector<double, 3UL> &v)
{
    return blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>{
        {0.00, -v[2UL], v[1UL]}, {v[2UL], 0.00, -v[0UL]}, {-v[1UL], v[0UL], 0.00}};
}

/**
 * @brief Computes the derivative of the quaternion orientation with respect to arc-length.
 *
 * @param u The local curvature vector of the backbone [1/m].
 * @param h The orientation quaternion of the local body frame (any dense vector expression).
 * @return The rate of change of the quaternion with respect to arc-length.
 */
template <typename VT>
inline blaze::StaticVector<double, 4UL> quaternionDiff(const blaze::StaticVector<double, 3UL> &u,
                                                       const blaze::DenseVector<VT, blaze::columnVector> &h_)
{
    const auto &h = *h_;
    return blaze::StaticVector<double, 4UL>{0.50 * (-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]),
                                            0.50 * (u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]),
                                            0.50 * (u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]),
                                            0.50 * (u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL])};
}

} // namespace ctr::math
