#pragma once

#include <blaze/Math.h>
#include <cmath>

namespace ctr
{

// ─── Compile-time robot constants ─────────────────────────────────────────────

inline constexpr std::size_t NUM_TUBES = 3UL;    ///< Number of concentric tubes.
inline constexpr std::size_t MAX_SEGMENTS = 7UL; ///< Max arc-length segments (2*NUM_TUBES + 1).
inline constexpr std::size_t BVP_DIM = 5UL;      ///< Dimension of the BVP initial guess / residue.

// ─── Type aliases ─────────────────────────────────────────────────────────────

/// ODE state vector: [mb_x, mb_y, uz_1, uz_2, uz_3, theta_1..3, pos_x..z, quat_w..z]
using state_type = blaze::StaticVector<double, 15UL>;

/// BVP boundary condition guess / residue vector.
using bvp_type = blaze::StaticVector<double, BVP_DIM>;

// ─── ODE state-vector index constants ─────────────────────────────────────────

namespace StateIdx
{
inline constexpr std::size_t MB_X = 0UL;    ///< Bending moment x (innermost tube, proximal)
inline constexpr std::size_t MB_Y = 1UL;    ///< Bending moment y (innermost tube, proximal)
inline constexpr std::size_t UZ_1 = 2UL;    ///< Torsional curvature, tube 1
inline constexpr std::size_t UZ_2 = 3UL;    ///< Torsional curvature, tube 2
inline constexpr std::size_t UZ_3 = 4UL;    ///< Torsional curvature, tube 3
inline constexpr std::size_t THETA_1 = 5UL; ///< Twist angle, tube 1
inline constexpr std::size_t THETA_2 = 6UL; ///< Twist angle, tube 2
inline constexpr std::size_t THETA_3 = 7UL; ///< Twist angle, tube 3
inline constexpr std::size_t POS_X = 8UL;   ///< Backbone position x
inline constexpr std::size_t POS_Y = 9UL;   ///< Backbone position y
inline constexpr std::size_t POS_Z = 10UL;  ///< Backbone position z
inline constexpr std::size_t QUAT_W = 11UL; ///< Orientation quaternion – scalar part
inline constexpr std::size_t QUAT_X = 12UL; ///< Orientation quaternion – x
inline constexpr std::size_t QUAT_Y = 13UL; ///< Orientation quaternion – y
inline constexpr std::size_t QUAT_Z = 14UL; ///< Orientation quaternion – z
} // namespace StateIdx

// ─── Shared BVP sanitization utility ──────────────────────────────────────────

namespace detail
{

/// Maximum magnitude of torsional curvature used to clamp an initial BVP guess.
inline constexpr double kMaxBVPTwist = 50.0;

/**
 * @brief Clamps a BVP initial guess to a physically reasonable range.
 *
 * Zeroes the moment components (indices 0-1) and clamps torsional-curvature
 * components (indices 2-4) to [-kMaxBVPTwist, kMaxBVPTwist].  Non-finite
 * values are replaced by 0.
 *
 * @param x BVP guess vector to sanitize in place.
 */
inline void sanitizeBVPGuess(bvp_type &x) noexcept
{
    for (std::size_t i = 2UL; i < BVP_DIM; ++i)
    {
        auto v = x[i];
        if (!std::isfinite(v))
            v = 0.0;
        x[i] = std::clamp(v, -kMaxBVPTwist, kMaxBVPTwist);
    }
    x[0UL] = 0.0;
    x[1UL] = 0.0;
}

} // namespace detail

} // namespace ctr
