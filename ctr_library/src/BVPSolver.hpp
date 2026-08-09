#pragma once

// PRIVATE header — the BVP solver strategy hierarchy is an implementation
// detail of the library and is not installed. Users select a solver through
// ctr::RootFindingMethod on the public CTR API.

#include "ctr/Types.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

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

/**
 * @brief Solves the small dense system A·x = b by Gaussian elimination with
 *        partial pivoting.
 *
 * Hand-rolled on purpose: blaze::solve/blaze::invert dispatch on the matrix
 * size at RUNTIME, so their generic LAPACK branches stay referenced in
 * unoptimized builds even though they are never taken for these fixed sizes —
 * which would silently reintroduce a LAPACK link dependency. For 3×3 and 5×5
 * systems this elimination is a handful of flops.
 *
 * @return The solution, or a NaN-filled vector if A is exactly singular.
 */
template <std::size_t N>
[[nodiscard]] inline blaze::StaticVector<double, N> solveLinear(Mat<N, N> A, blaze::StaticVector<double, N> b) noexcept
{
    for (std::size_t k = 0UL; k < N; ++k)
    {
        // Partial pivoting.
        std::size_t p = k;
        for (std::size_t i = k + 1UL; i < N; ++i)
            if (std::fabs(A(i, k)) > std::fabs(A(p, k)))
                p = i;
        if (p != k)
        {
            for (std::size_t j = k; j < N; ++j)
                std::swap(A(k, j), A(p, j));
            std::swap(b[k], b[p]);
        }

        const double piv = A(k, k);
        if (piv == 0.0)
            return blaze::StaticVector<double, N>(std::numeric_limits<double>::quiet_NaN());

        for (std::size_t i = k + 1UL; i < N; ++i)
        {
            const double m = A(i, k) / piv;
            if (m != 0.0)
            {
                for (std::size_t j = k + 1UL; j < N; ++j)
                    A(i, j) -= m * A(k, j);
                b[i] -= m * b[k];
            }
        }
    }

    // Back substitution.
    blaze::StaticVector<double, N> x;
    for (std::size_t k = N; k-- > 0UL;)
    {
        double s = b[k];
        for (std::size_t j = k + 1UL; j < N; ++j)
            s -= A(k, j) * x[j];
        x[k] = s / A(k, k);
    }
    return x;
}

/// Inverts a small dense matrix column-by-column via solveLinear.
template <std::size_t N> [[nodiscard]] inline Mat<N, N> invertMatrix(const Mat<N, N> &A) noexcept
{
    Mat<N, N> inv;
    blaze::StaticVector<double, N> e(0.0);
    for (std::size_t j = 0UL; j < N; ++j)
    {
        e[j] = 1.0;
        blaze::column(inv, j) = solveLinear<N>(A, e);
        e[j] = 0.0;
    }
    return inv;
}

/**
 * @brief Inverts a BVP Jacobian.
 *
 * Direct pivoted inversion (allocation- and LAPACK-free). Near singularity it
 * falls back to the damped pseudo-inverse (JᵀJ + μI)⁻¹Jᵀ with μ scaled to the
 * Jacobian's magnitude; if even that is non-finite, returns a zero matrix
 * (turning the caller's step into a detectable stall rather than a NaN
 * cascade).
 */
[[nodiscard]] inline Mat<BVP_DIM, BVP_DIM> invertJacobian(const Mat<BVP_DIM, BVP_DIM> &J)
{
    const Mat<BVP_DIM, BVP_DIM> Jinv = invertMatrix<BVP_DIM>(J);
    if (blaze::isfinite(Jinv))
        return Jinv;

    Mat<BVP_DIM, BVP_DIM> A = blaze::trans(J) * J;
    double mu = 1.0e-10 * blaze::max(blaze::abs(blaze::diagonal(A)));
    if (!(mu > 0.0) || !std::isfinite(mu))
        mu = 1.0e-10;
    for (std::size_t i = 0UL; i < BVP_DIM; ++i)
        A(i, i) += mu;

    const Mat<BVP_DIM, BVP_DIM> pinv = invertMatrix<BVP_DIM>(A) * blaze::trans(J);
    if (blaze::isfinite(pinv))
        return pinv;
    return Mat<BVP_DIM, BVP_DIM>(0.0);
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
