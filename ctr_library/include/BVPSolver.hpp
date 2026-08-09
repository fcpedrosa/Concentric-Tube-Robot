#pragma once

#include "CTRTypes.hpp"
#include "mathOperations.hpp"
#include <memory>

namespace ctr {

class CTR; // forward declaration — avoids circular #include with CTR.hpp

// ─── Abstract base ─────────────────────────────────────────────────────────────

/**
 * @brief Strategy interface for BVP root-finding algorithms.
 *
 * Each concrete solver implements the shooting-method loop that drives the
 * BVP residue (returned by CTR::ODESolver) to zero.
 */
class BVPSolver
{
public:
    virtual ~BVPSolver() = default;

    /**
     * @brief Solves the CTR boundary value problem.
     *
     * @param initGuess Proximal boundary condition vector (updated in place on convergence).
     * @param ctr       CTR object that provides ODE residue and Jacobian evaluations.
     * @return true if the BVP converged within the iteration limit, false otherwise.
     */
    [[nodiscard]] virtual bool solve(bvp_type &initGuess, CTR &ctr) = 0;
};

// ─── Concrete solvers ──────────────────────────────────────────────────────────

/// Powell Dog-Leg trust-region solver.
class PowellDogLegSolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

/// Levenberg-Marquardt damped-least-squares solver.
class LevenbergMarquardtSolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

/// Broyden rank-1 quasi-Newton solver (inverse Jacobian update).
class BroydenSolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

/// Broyden rank-1 quasi-Newton solver (direct Jacobian update).
class BroydenIISolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

/// Classic Newton-Raphson solver with PID-damped update.
class NewtonRaphsonSolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

/// Modified Newton-Raphson solver with Armijo line search.
class ModifiedNewtonRaphsonSolver final : public BVPSolver
{
public:
    [[nodiscard]] bool solve(bvp_type &initGuess, CTR &ctr) override;
};

// ─── Factory ──────────────────────────────────────────────────────────────────

/**
 * @brief Creates a BVPSolver of the requested type.
 *
 * @param method Root-finding algorithm to use.
 * @return Owning pointer to a freshly created concrete BVPSolver.
 */
[[nodiscard]] std::unique_ptr<BVPSolver> makeBVPSolver(mathOp::rootFindingMethod method);

} // namespace ctr
