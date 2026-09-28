#pragma once

/**
 * @file Types.hpp
 * @brief Library-wide constants, vector/matrix aliases, solver selection and
 *        the result/option types returned by ctr::CTR.
 */

#include <blaze/Math.h>
#include <cstdint>

/**
 * @brief Everything in the CTR kinematics library.
 *
 * Start with ctr::CTR (the robot) and ctr::Tube (its components). All
 * vectors and matrices are Blaze fixed-size types; see @ref conventions for
 * units, frames and vector layouts.
 */
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

/// ODE state vector: [mb_x, mb_y, uz_1, uz_2, uz_3, theta_1..3, pos_x..z, quat_w..z].
/// Index it with the StateIdx constants; see @ref conv_state for units.
using state_type = blaze::StaticVector<double, 15UL>;

/// Shooting vector of the forward-kinematics BVP (the five unknown proximal
/// values, non-dimensionalized to [1/m]; see @ref conv_bvp).
///
/// Users only need to create one (`ctr::bvp_type guess{};` — zero is a valid
/// cold start) and keep passing it to CTR::actuate / CTR::solveIK, which
/// update it in place for warm starting.
using bvp_type = blaze::StaticVector<double, BVP_DIM>;

/// Fixed-size M×N matrix of doubles in the library-wide storage order
/// (column-major), e.g. the 3×6 `Mat<3, 6>` returned by CTR::kinematicJacobian.
/// Use this alias for every fixed-size matrix so no storage-order-converting
/// copies sneak in at API boundaries.
template <std::size_t M, std::size_t N> using Mat = blaze::StaticMatrix<double, M, N, blaze::columnMajor>;

// ─── ODE state-vector index constants ─────────────────────────────────────────

/// Indices into ctr::state_type (see @ref conv_state for units and frames).
namespace StateIdx
{
inline constexpr std::size_t MB_X = 0UL;    ///< Transverse bending moment x [N·m], tube-1 body frame
inline constexpr std::size_t MB_Y = 1UL;    ///< Transverse bending moment y [N·m], tube-1 body frame
inline constexpr std::size_t UZ_1 = 2UL;    ///< Torsional curvature, tube 1 [1/m]
inline constexpr std::size_t UZ_2 = 3UL;    ///< Torsional curvature, tube 2 [1/m]
inline constexpr std::size_t UZ_3 = 4UL;    ///< Torsional curvature, tube 3 [1/m]
inline constexpr std::size_t THETA_1 = 5UL; ///< Twist angle, tube 1 [rad] (≡ 0: the reference)
inline constexpr std::size_t THETA_2 = 6UL; ///< Twist angle of tube 2 relative to tube 1 [rad]
inline constexpr std::size_t THETA_3 = 7UL; ///< Twist angle of tube 3 relative to tube 1 [rad]
inline constexpr std::size_t POS_X = 8UL;   ///< Backbone position x [m], global frame
inline constexpr std::size_t POS_Y = 9UL;   ///< Backbone position y [m], global frame
inline constexpr std::size_t POS_Z = 10UL;  ///< Backbone position z [m], global frame
inline constexpr std::size_t QUAT_W = 11UL; ///< Orientation quaternion – scalar part
inline constexpr std::size_t QUAT_X = 12UL; ///< Orientation quaternion – x
inline constexpr std::size_t QUAT_Y = 13UL; ///< Orientation quaternion – y
inline constexpr std::size_t QUAT_Z = 14UL; ///< Orientation quaternion – z
} // namespace StateIdx

// ─── Solver selection ─────────────────────────────────────────────────────────

/**
 * @brief Root-finding algorithm used by the shooting method to solve the
 *        CTR's boundary value problem.
 *
 * The default is the right choice for almost every application: it is the
 * fastest on warm starts and falls back automatically to the more robust
 * methods on hard cold starts. The others are exposed mainly for
 * cross-validation and experimentation. Select one in the CTR constructor or
 * with CTR::setBVPMethod.
 */
enum class RootFindingMethod : std::uint8_t
{
    ModifiedNewtonRaphson, ///< Newton iteration with Armijo line search (Stoer & Bulirsch); falls back to
                           ///< Broyden, dog-leg, then Levenberg-Marquardt on hard cold starts. Default.
    LevenbergMarquardt,    ///< Damped least squares with adaptive damping.
    PowellDogLeg,          ///< Trust-region dog-leg method.
    Broyden                ///< Rank-1 quasi-Newton (inverse-Jacobian update).
};

// ─── Result types ─────────────────────────────────────────────────────────────

/// Termination status of a BVP solve (see FKResult::status, IKResult::lastBVPStatus).
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
 *
 * @code{.cpp}
 * if (const ctr::FKResult fk = robot.actuate(q, guess); !fk)
 *     std::cerr << "FK failed after " << fk.iterations << " iterations\n";
 * @endcode
 */
struct FKResult
{
    SolverStatus status{SolverStatus::NumericalError}; ///< How the solver terminated.
    std::size_t iterations{0UL};                       ///< Solver iterations used.
    double residual{0.0};                              ///< Final non-dimensional BVP residue norm (L∞) [1/m].

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
 *
 * The step caps make the iteration count scale with joint-space DISTANCE, not
 * just with the local convergence rate: a target that needs Δβ of translation
 * costs at least max|Δβ| / maxBetaStep iterations (and max|Δα| / maxAlphaStep
 * for rotation) before Levenberg-Marquardt damping shrinks steps below the
 * caps. maxIterations must therefore cover the excursion, not merely the
 * endgame — the default accommodates a full-stroke β retraction/extension.
 * Lower it for latency-bounded tracking loops, where a warm-started solve over
 * a millimetre-scale target step converges in well under ten iterations.
 *
 * @code{.cpp}
 * ctr::IKOptions opts;
 * opts.maxIterations = 20;          // real-time tracking: bounded latency
 * auto ik = robot.solveIK(target, 5e-4, guess, opts);
 * @endcode
 */
struct IKOptions
{
    std::size_t maxIterations = 200UL; ///< Iteration budget (Jacobian evaluations).
    std::size_t maxBacktracks = 8UL;   ///< Step halvings tried per iteration before re-damping.
    double maxBetaStep = 2.0e-3;       ///< Per-iteration cap on each |Δβ| [m].
    double maxAlphaStep = 0.35;        ///< Per-iteration cap on each |Δα| [rad].
    double dampingSeed = 1.0e-3;       ///< LM damping floor/seed, relative to max(diag(JJᵀ)).
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
    bool converged{false};                                    ///< ||tip - target|| <= posTol at exit.
    double positionError{0.0};                                ///< Final ||tip - target|| [m].
    std::size_t iterations{0UL};                              ///< IK iterations used.
    SolverStatus lastBVPStatus{SolverStatus::NumericalError}; ///< Status of the last internal BVP solve.
    blaze::StaticVector<double, 6UL> q{};                     ///< Joint configuration at exit.

    /** @brief True iff the tip reached the target within posTol. */
    [[nodiscard]] explicit operator bool() const noexcept { return converged; }
};

} // namespace ctr
