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
 * @brief Computes the closest congruent of an angle within the [0, 3Pi] interval.
 *
 * @param angle Angle in radians.
 * @return The corresponding closest congruent angle in radians within [0, 3Pi].
 */
inline double congruentAngle(double angle)
{
    constexpr double THREE_PI = 3.00 * std::numbers::pi;
    return std::fmod(std::fabs(angle), THREE_PI) * (angle < 0.00 ? -1.00 : 1.00);
}

/**
 * @brief Computes the damped least-squares pseudo-inverse of a DynamicMatrix via SVD.
 *
 * Heavy implementation lives in mathOperations.cpp to avoid bloating every
 * translation unit that includes this header.
 *
 * @param M Input column-major DynamicMatrix.
 * @return  The pseudo-inverse as a column-major DynamicMatrix.
 */
blaze::DynamicMatrix<double, blaze::columnMajor> pInv(blaze::DynamicMatrix<double, blaze::columnMajor> M);

/**
 * @brief Thin wrapper: accepts any Blaze matrix type and delegates to the
 *        non-template overload above.
 *
 * @tparam MT  Blaze matrix type (deduced).
 * @tparam SO  Storage order (deduced).
 * @param  M   The input matrix.
 * @return     The pseudo-inverse as a column-major DynamicMatrix.
 */
template <typename MT, bool SO> inline auto pInv(const blaze::Matrix<MT, SO> &M)
{
    return pInv(blaze::DynamicMatrix<double, blaze::columnMajor>(*M));
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
 * @brief Computes the pre-multiplication of a 3×3 matrix by the hat operator of a
 * 3-dimensional vector, i.e. the operation \hat{v} * M, without forming \hat{v}.
 *
 * Accepts any dense-matrix expression (e.g. blaze::trans of a rotation matrix)
 * so callers do not materialize a temporary.
 *
 * @param v The input 3-dimensional vector.
 * @param M The 3×3 matrix expression to be pre-multiplied.
 * @return The 3×3 matrix result of the operation.
 */
template <typename MT, bool SO>
inline blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>
hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::Matrix<MT, SO> &M_)
{
    const auto &M = *M_;
    return blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>{
        {-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL],
         -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
        {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL],
         M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
        {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL],
         -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};
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
