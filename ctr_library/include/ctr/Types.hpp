#pragma once

#include <blaze/Math.h>
#include <cstdint>

namespace ctr
{

// ─── Compile-time robot constants ─────────────────────────────────────────────

inline constexpr std::size_t NUM_TUBES = 3UL; ///< Number of concentric tubes.

/// Maximum number of tube-transition points along the backbone (each of the
/// NUM_TUBES tubes contributes a curved-section start and a distal end, plus
/// the robot base at s = 0). The number of arc-length *segments* is at most
/// MAX_SEGMENTS - 1.
inline constexpr std::size_t MAX_SEGMENTS = 7UL;

inline constexpr std::size_t BVP_DIM = 5UL; ///< Dimension of the BVP initial guess / residue.

// ─── Type aliases ─────────────────────────────────────────────────────────────

/// ODE state vector: [mb_x, mb_y, uz_1, uz_2, uz_3, theta_1..3, pos_x..z, quat_w..z]
using state_type = blaze::StaticVector<double, 15UL>;

/// BVP boundary condition guess / residue vector.
using bvp_type = blaze::StaticVector<double, BVP_DIM>;

/// Fixed-size matrix in the library-wide storage order (column-major).
/// Use this alias for every fixed-size matrix so no storage-order-converting
/// copies sneak in at API boundaries.
template <std::size_t M, std::size_t N> using Mat = blaze::StaticMatrix<double, M, N, blaze::columnMajor>;

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

// ─── Solver selection ─────────────────────────────────────────────────────────

/**
 * @brief Root-finding algorithm used by the shooting method to solve the
 *        CTR's boundary value problem.
 */
enum class RootFindingMethod : std::uint8_t
{
    ModifiedNewtonRaphson, ///< Newton iteration with Armijo line search (Stoer & Bulirsch). Default.
    LevenbergMarquardt,    ///< Damped least squares with adaptive damping.
    PowellDogLeg,          ///< Trust-region dog-leg method.
    Broyden                ///< Rank-1 quasi-Newton (inverse-Jacobian update).
};

// ─── Result types ─────────────────────────────────────────────────────────────

/// Termination status of a BVP solve.
enum class SolverStatus : std::uint8_t
{
    Converged,     ///< Residue norm fell below the tolerance.
    MaxIterations, ///< Iteration budget exhausted before convergence.
    NumericalError ///< NaN/Inf or a singular update was encountered.
};

/**
 * @brief Outcome of a forward-kinematics (BVP) solve.
 *
 * Contextually convertible to bool: true iff the BVP converged.
 */
struct FKResult
{
    SolverStatus status{SolverStatus::NumericalError}; ///< How the solver terminated.
    std::size_t iterations{0UL};                       ///< Solver iterations used.
    double residual{0.0};                              ///< Final BVP residue norm (L∞).

    /** @brief True iff the BVP converged. */
    [[nodiscard]] explicit operator bool() const noexcept { return status == SolverStatus::Converged; }
};

/**
 * @brief Tuning knobs for the inverse-kinematics solver.
 *
 * The defaults are appropriate for tabletop-scale CTRs (backbone lengths of
 * tens of centimeters); they bound the per-iteration step so the damped
 * least-squares iteration stays inside the region where its linear model and
 * the warm-started BVP solve are reliable.
 */
struct IKOptions
{
    std::size_t maxIterations = 50UL; ///< Iteration budget.
    double maxBetaStep = 2.0e-3;      ///< Per-iteration cap on each |Δβ| [m].
    double maxAlphaStep = 0.35;       ///< Per-iteration cap on each |Δα| [rad].
    double dampingSeed = 1.0e-3;      ///< Initial LM damping relative to max(diag(JJᵀ)).
};

/**
 * @brief Outcome of an inverse-kinematics solve.
 *
 * Contextually convertible to bool: true iff the tip reached the target
 * within the requested position tolerance (NOT merely whether the last
 * internal BVP solve converged).
 */
struct IKResult
{
    bool converged{false};                                   ///< ||tip - target|| <= posTol at exit.
    double positionError{0.0};                               ///< Final ||tip - target|| [m].
    std::size_t iterations{0UL};                             ///< IK iterations used.
    SolverStatus lastBVPStatus{SolverStatus::NumericalError}; ///< Status of the last internal BVP solve.
    blaze::StaticVector<double, 6UL> q{};                    ///< Joint configuration at exit.

    /** @brief True iff the tip reached the target within posTol. */
    [[nodiscard]] explicit operator bool() const noexcept { return converged; }
};

} // namespace ctr
