#include "BVPSolver.hpp"
#include "ctr/CTR.hpp"

#include <cmath>
#include <limits>

namespace ctr
{

// ─── Utility shared by all solvers ───────────────────────────────────────────

using detail::sanitizeBVPGuess;

namespace
{

/// Assembles the FKResult for a finished solve.
[[nodiscard]] FKResult makeResult(bool found, std::size_t iterations, const bvp_type &residue) noexcept
{
    const double res_norm = blaze::linfNorm(residue);
    SolverStatus status = SolverStatus::MaxIterations;
    if (found)
        status = SolverStatus::Converged;
    else if (!blaze::isfinite(residue))
        status = SolverStatus::NumericalError;
    return {status, iterations, res_norm};
}

/// Shared iteration budget: warm-started solves need < 10 iterations; a solve
/// that has not converged by 50 will not converge (caps worst-case cost).
inline constexpr std::size_t kMaxSolverIterations = 50UL;

} // namespace

// ─── Powell Dog-Leg ──────────────────────────────────────────────────────────

FKResult PowellDogLegSolver::solve(bvp_type &initGuess, ShootingProblem &problem)
{
    bool found;
    std::size_t k = 0UL;
    constexpr std::size_t k_max = kMaxSolverIterations;
    double alpha, beta, delta, eps1, eps2, rho, c;
    bvp_type g, f, f_new, x_new, h_sd, h_gn, h_dl;
    Mat<BVP_DIM, BVP_DIM> J;

    sanitizeBVPGuess(initGuess);

    delta = 1.0;
    eps1 = eps2 = 1.0e-12;

    f = problem.residual(initGuess);
    J = problem.jacobian(initGuess, f);
    g = blaze::trans(J) * f;

    found = ((blaze::linfNorm(f) <= problem.tolerance()) || (blaze::linfNorm(g) <= eps1));

    while (!found && (k < k_max))
    {
        k++;

        alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
        h_sd = -alpha * g;
        h_gn = -(detail::invertJacobian(J) * f);

        if (blaze::norm(h_gn) <= delta)
            h_dl = h_gn;
        else
        {
            if (blaze::norm(h_sd) >= delta)
                h_dl = delta * blaze::normalize(h_sd);
            else
            {
                c = blaze::trans(h_sd) * (h_gn - h_sd);
                if (c <= 0.0)
                    beta =
                        (-c + std::sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) /
                        blaze::sqrNorm(h_gn - h_sd);
                else
                    beta =
                        (delta * delta - blaze::sqrNorm(h_sd)) /
                        (c + std::sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
                h_dl = h_sd + beta * (h_gn - h_sd);
            }
        }

        if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
            found = true;
        else
        {
            x_new = initGuess + h_dl;
            f_new = problem.residual(x_new);
            // Gain ratio = actual reduction / model-predicted reduction, with the
            // Gauss-Newton model L(h) = ||f + J h||²: predicted = ||f||² − ||f + J h_dl||².
            const double predictedReduction = blaze::sqrNorm(f) - blaze::sqrNorm(f + J * h_dl);
            rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) /
                  std::max(predictedReduction, std::numeric_limits<double>::min());

            if (rho > 0.0)
            {
                initGuess = std::move(x_new);
                f = std::move(f_new);
                J = problem.jacobian(initGuess, f);
                g = blaze::trans(J) * f;
                if ((blaze::linfNorm(f) <= problem.tolerance()) || (blaze::linfNorm(g) <= eps1))
                    found = true;
            }

            if (rho > 0.75)
                delta = std::max(delta, 3.0 * blaze::norm(h_dl));
            else if (rho < 0.25)
                delta *= 0.5;

            if (delta < eps2 * (blaze::norm(initGuess) + eps2))
                found = true;
        }
    }

    return makeResult(found, k, f);
}

// ─── Levenberg-Marquardt ─────────────────────────────────────────────────────

FKResult LevenbergMarquardtSolver::solve(bvp_type &initGuess, ShootingProblem &problem)
{
    std::size_t k = 0UL;
    constexpr std::size_t k_max = kMaxSolverIterations;
    bvp_type h, g, f, f_new;
    Mat<BVP_DIM, BVP_DIM> J, A;
    blaze::IdentityMatrix<double> I(BVP_DIM);
    double rho, nu = 2.0, mu, tau = 1.0e-3, e1 = 1.0e-12;
    bool found;

    sanitizeBVPGuess(initGuess);

    f = problem.residual(initGuess);
    J = problem.jacobian(initGuess, f);
    A = blaze::trans(J) * J;
    g = blaze::trans(J) * f;
    found = (blaze::linfNorm(g) <= e1);
    mu = tau * blaze::max(blaze::diagonal(A));

    while (!found && (k < k_max))
    {
        k++;
        h = detail::solveLinear<BVP_DIM>(A + (mu * I), -g);

        f_new = problem.residual(initGuess + h);
        rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h) * ((mu * h) - g));

        if (rho > 0.0)
        {
            initGuess += h;
            J = problem.jacobian(initGuess, f_new);
            A = blaze::trans(J) * J;
            f = std::move(f_new);
            g = blaze::trans(J) * f;
            found = (blaze::linfNorm(g) <= e1);
            mu = mu * std::max(0.33333333, 1.0 - blaze::pow(2.0 * rho - 1.0, 3.0));
            nu = 2.0;
        }
        else
        {
            mu = mu * nu;
            nu = 2.0 * nu;
        }

        if (blaze::linfNorm(f) <= problem.tolerance())
            found = true;
    }

    return makeResult(found, k, f);
}

// ─── Broyden ─────────────────────────────────────────────────────────────────

FKResult BroydenSolver::solve(bvp_type &initGuess, ShootingProblem &problem)
{
    sanitizeBVPGuess(initGuess);

    Mat<BVP_DIM, BVP_DIM> JacInv;
    bvp_type F, Fold, X, Xold, deltaX, deltaF;

    X = initGuess;
    F = problem.residual(X);
    JacInv = detail::invertJacobian(problem.jacobian(X, F));
    bool found = (blaze::linfNorm(F) <= problem.tolerance());

    std::size_t k = 0UL;
    constexpr std::size_t k_max = kMaxSolverIterations;
    while (!found && (k < k_max))
    {
        k++;

        // Quasi-Newton step with the CURRENT inverse Jacobian.
        Xold = X;
        Fold = F;
        X = Xold - JacInv * Fold;
        F = problem.residual(X);

        // Bounded recovery from a non-finite residual.
        std::size_t recovery = 0UL;
        while (!blaze::isfinite(F))
        {
            if (++recovery > 20UL)
            {
                initGuess = std::move(X);
                return {SolverStatus::NumericalError, k, blaze::linfNorm(Fold)};
            }
            X *= 0.75;
            sanitizeBVPGuess(X);
            F = problem.residual(X);
            JacInv = detail::invertJacobian(problem.jacobian(X, F));
        }

        if (k % 10UL == 0UL)
        {
            // Periodic exact refresh keeps the update from drifting.
            JacInv = detail::invertJacobian(problem.jacobian(X, F));
        }
        else
        {
            // Sherman-Morrison ("good Broyden") inverse update:
            //   JacInv += (dx - JacInv*dF) dxᵀ JacInv / (dxᵀ JacInv dF)
            deltaX = X - Xold;
            deltaF = F - Fold;
            const double denom = blaze::trans(deltaX) * (JacInv * deltaF);
            if ((blaze::sqrNorm(deltaX) > 0.0) && (std::fabs(denom) > 1.0e-300))
                JacInv += ((deltaX - JacInv * deltaF) / denom) * (blaze::trans(deltaX) * JacInv);
        }

        found = (blaze::linfNorm(F) <= problem.tolerance());
    }

    initGuess = std::move(X);
    return makeResult(found, k, F);
}

// ─── Modified Newton-Raphson ─────────────────────────────────────────────────

FKResult ModifiedNewtonRaphsonSolver::solve(bvp_type &initGuess, ShootingProblem &problem)
{
    /*
     * Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition
     * by Josef Stoer & Roland Bulirsch.
     */
    sanitizeBVPGuess(initGuess);

    bool found;
    bvp_type f(problem.residual(initGuess)), d;
    blaze::StaticVector<double, BVP_DIM, blaze::rowVector> Dh;
    Mat<BVP_DIM, BVP_DIM> D, D_inv;
    double h, h_0, gamma, improvementFactor, d_norm, Dh_norm;
    std::size_t j = 0UL, k = 0UL;
    constexpr std::size_t k_max = kMaxSolverIterations;
    constexpr std::size_t j_max = 60UL; // 0.5^60 ~ 1e-18: step is numerically nil

    found = (blaze::linfNorm(f) <= problem.tolerance());

    auto setupMethod = [&]() -> void
    {
        f = problem.residual(initGuess);
        D = problem.jacobian(initGuess, f);

        std::size_t isfinite_guard = 0UL;
        while (!blaze::isfinite(D))
        {
            initGuess *= 0.75;
            sanitizeBVPGuess(initGuess);
            f = problem.residual(initGuess);
            D = problem.jacobian(initGuess, f);
            if (++isfinite_guard > 20UL)
                break;
        }

        D_inv = detail::invertJacobian(D);
        d = D_inv * f;
        gamma = 1.0 / (blaze::norm(D_inv) * blaze::norm(D));
        h_0 = blaze::sqrNorm(f);
        Dh = 2.0 * blaze::trans(f) * D;
        d_norm = blaze::norm(d);
        Dh_norm = blaze::norm(Dh);
    };

    while (!found && (k < k_max))
    {
        k++;
        setupMethod();
        j = 0UL;

        // Armijo backtracking: find the largest step 0.5^j that sufficiently
        // decreases ||f||². The accepted step below uses exactly this j, so the
        // residue f (and the recorded backbone trajectory) always corresponds
        // to the point the iterate actually moves to.
        while (true)
        {
            f = problem.residual(initGuess - blaze::pow(0.5, j) * d);

            while (!blaze::isfinite(f))
            {
                j++;
                f = problem.residual(initGuess - blaze::pow(0.5, j) * d);
                if (j > 10UL)
                {
                    initGuess *= 0.75;
                    sanitizeBVPGuess(initGuess);
                    setupMethod();
                    j = 0UL;
                }
            }

            h = blaze::sqrNorm(f);
            improvementFactor = blaze::pow(0.5, j) * 0.25 * gamma * d_norm * Dh_norm;

            if ((h <= (h_0 - improvementFactor)) || (j >= j_max))
                break;
            j++;
        }

        initGuess -= blaze::pow(0.5, j) * d;

        if (blaze::linfNorm(f) <= problem.tolerance())
            return makeResult(true, k, f);
    }

    if (!found)
    {
        // Fallback cascade. Broyden goes first: as a quasi-Newton method it
        // follows a different trajectory than the Newton-type solvers and is
        // empirically able to escape the non-root local minima of ||f||² that
        // trap MNR/PDL/LM on rotation-heavy cold starts. Each stage restarts
        // from a progressively cleaner guess — continuing from the stalled
        // iterate keeps a solver trapped in the same basin, so the later
        // stages restart from zero outright.
        sanitizeBVPGuess(initGuess);
        initGuess *= 0.75;
        FKResult fallback = BroydenSolver{}.solve(initGuess, problem);

        if (!fallback)
        {
            initGuess = 0.0; // clean cold restart, away from the stalled basin
            fallback = BroydenSolver{}.solve(initGuess, problem);
        }

        if (!fallback)
        {
            initGuess = 0.0;
            fallback = PowellDogLegSolver{}.solve(initGuess, problem);
        }

        if (!fallback)
        {
            initGuess = 0.0;
            fallback = LevenbergMarquardtSolver{}.solve(initGuess, problem);
        }

        fallback.iterations += k;
        return fallback;
    }

    return makeResult(found, k, f);
}

// ─── Factory ─────────────────────────────────────────────────────────────────

std::unique_ptr<BVPSolver> makeBVPSolver(RootFindingMethod method)
{
    switch (method)
    {
    case RootFindingMethod::LevenbergMarquardt:
        return std::make_unique<LevenbergMarquardtSolver>();
    case RootFindingMethod::PowellDogLeg:
        return std::make_unique<PowellDogLegSolver>();
    case RootFindingMethod::ModifiedNewtonRaphson:
        return std::make_unique<ModifiedNewtonRaphsonSolver>();
    case RootFindingMethod::Broyden:
        return std::make_unique<BroydenSolver>();
    }
    return std::make_unique<ModifiedNewtonRaphsonSolver>(); // unreachable, appeases the compiler
}

} // namespace ctr
