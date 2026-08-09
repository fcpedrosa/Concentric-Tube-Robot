#pragma once

// PRIVATE header — the BVP solver strategy hierarchy is an implementation
// detail of the library and is not installed. Users select a solver through
// ctr::RootFindingMethod on the public CTR API.

#include "ctr/Types.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

namespace ctr
{

class CTR; // forward declaration — avoids circular #include with CTR.hpp

// ─── Shared BVP-guess sanitization ─────────────────────────────────────────────

namespace detail
{

/// Maximum magnitude of any component of a BVP guess. The shooting vector is
/// non-dimensionalized (all components are curvatures [1/m]), so one bound
/// applies uniformly.
inline constexpr double kMaxBVPTwist = 50.0;

/**
 * @brief Clamps a BVP initial guess to a physically reasonable range.
 *
 * Non-finite values are replaced by 0; every component is clamped to
 * ±kMaxBVPTwist (all components are curvatures [1/m]).
 *
 * The moment components are deliberately NOT zeroed: with a distal load the
 * proximal bending moment is nonzero, and zeroing would both make loaded
 * BVPs unsolvable and destroy warm starts.
 */
inline void sanitizeBVPGuess(bvp_type &x) noexcept
{
    for (std::size_t i = 0UL; i < BVP_DIM; ++i)
    {
        auto v = x[i];
        if (!std::isfinite(v))
            v = 0.0;
        x[i] = std::clamp(v, -kMaxBVPTwist, kMaxBVPTwist);
    }
}

} // namespace detail

// ─── Shooting-problem facade ───────────────────────────────────────────────────

/**
 * @brief Narrow interface through which BVP solvers evaluate the shooting problem.
 *
 * Constructible only by CTR (the sole friend). Keeps the FK internals
 * (residual evaluation, Jacobians, integration buffers) private on CTR while
 * still letting the solver strategies drive them.
 */
class ShootingProblem
{
    friend class CTR;

    explicit ShootingProblem(CTR &ctr) noexcept : m_ctr(ctr) {}

    CTR &m_ctr;

  public:
    /** @brief One forward shot: integrates the ODE and returns the BVP residue. */
    [[nodiscard]] bvp_type residual(const bvp_type &x);

    /** @brief Finite-difference Jacobian of the residue w.r.t. the shooting variables. */
    [[nodiscard]] Mat<BVP_DIM, BVP_DIM> jacobian(const bvp_type &x, const bvp_type &f0);

    /** @brief BVP convergence tolerance (L∞ norm on the residue). */
    [[nodiscard]] double tolerance() const noexcept;
};

// ─── Abstract base ─────────────────────────────────────────────────────────────

/**
 * @brief Strategy interface for BVP root-finding algorithms.
 *
 * Each concrete solver implements the shooting-method loop that drives the
 * BVP residue to zero. Solvers are stateless.
 */
class BVPSolver
{
  public:
    virtual ~BVPSolver() = default;

    /**
     * @brief Solves the CTR boundary value problem.
     *
     * @param initGuess Proximal boundary condition vector (updated in place on convergence).
     * @param problem   Facade for residue/Jacobian evaluations.
     * @return FKResult describing convergence, iterations used, and final residue norm.
     */
    [[nodiscard]] virtual FKResult solve(bvp_type &initGuess, ShootingProblem &problem) = 0;
};

// ─── Concrete solvers ──────────────────────────────────────────────────────────

/// Powell Dog-Leg trust-region solver.
class PowellDogLegSolver final : public BVPSolver
{
  public:
    [[nodiscard]] FKResult solve(bvp_type &initGuess, ShootingProblem &problem) override;
};

/// Levenberg-Marquardt damped-least-squares solver.
class LevenbergMarquardtSolver final : public BVPSolver
{
  public:
    [[nodiscard]] FKResult solve(bvp_type &initGuess, ShootingProblem &problem) override;
};

/// Broyden rank-1 quasi-Newton solver (inverse Jacobian update).
class BroydenSolver final : public BVPSolver
{
  public:
    [[nodiscard]] FKResult solve(bvp_type &initGuess, ShootingProblem &problem) override;
};

/// Modified Newton-Raphson solver with Armijo line search (Stoer & Bulirsch).
class ModifiedNewtonRaphsonSolver final : public BVPSolver
{
  public:
    [[nodiscard]] FKResult solve(bvp_type &initGuess, ShootingProblem &problem) override;
};

// ─── Factory ──────────────────────────────────────────────────────────────────

/**
 * @brief Creates a BVPSolver of the requested type.
 *
 * @param method Root-finding algorithm to use.
 * @return Owning pointer to a freshly created concrete BVPSolver.
 */
[[nodiscard]] std::unique_ptr<BVPSolver> makeBVPSolver(RootFindingMethod method);

} // namespace ctr
