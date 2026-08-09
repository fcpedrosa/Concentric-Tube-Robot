#include "BVPSolver.hpp"
#include "ctr/CTR.hpp"
#include "ctr/detail/mathOperations.hpp"

#include <cmath>

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

} // namespace

// ─── Powell Dog-Leg ──────────────────────────────────────────────────────────

FKResult PowellDogLegSolver::solve(bvp_type &initGuess, ShootingProblem &problem)
{
    bool found;
    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    double alpha, beta, delta, eps1, eps2, rho, c;
    bvp_type g, f, f_new, x_new, h_sd, h_gn, h_dl;
    Mat<BVP_DIM, BVP_DIM> J;

    sanitizeBVPGuess(initGuess);

    delta = 1.0;
    eps1 = eps2 = 1.0e-22;

    f = problem.residual(initGuess);
    J = problem.jacobian(initGuess, f);
    g = blaze::trans(J) * f;

    found = ((blaze::linfNorm(f) <= problem.tolerance()) || (blaze::linfNorm(g) <= eps1));

    while (!found && (k < k_max))
    {
        k++;

        alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
        h_sd = -alpha * g;
        h_gn = -math::pInv(J) * f;

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
            rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h_dl) * ((delta * h_dl) - g));

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
    constexpr std::size_t k_max = 300UL;
    bvp_type h, g, f, f_new;
    Mat<BVP_DIM, BVP_DIM> J, A;
    blaze::IdentityMatrix<double> I(BVP_DIM);
    double rho, nu = 2.0, mu, tau = 1.0e-3, e1 = 1.0e-18;
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
        blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

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
    bool found;
    sanitizeBVPGuess(initGuess);

    Mat<BVP_DIM, BVP_DIM> JacInv, JacInvNew;
    bvp_type F, Fold, X, Xold, deltaX, deltaF;

    F = problem.residual(initGuess);
    X = std::move(initGuess);
    JacInvNew = JacInv = math::pInv(problem.jacobian(X, F));
    found = (blaze::linfNorm(F) <= problem.tolerance());

    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    while (!found && (k < k_max))
    {
        k++;

        deltaX = X - Xold;
        deltaF = F - Fold;

        JacInv = std::move(JacInvNew);
        if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.0))
            JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) *
                                     blaze::trans(deltaX) * JacInv;
        else
            JacInvNew = JacInv;

        Xold = std::move(X);
        Fold = std::move(F);
        X = Xold - JacInv * F;
        F = problem.residual(X);

        while (blaze::isnan(F))
        {
            X *= 0.75;
            sanitizeBVPGuess(X);
            F = problem.residual(X);
            JacInv = JacInvNew = math::pInv(problem.jacobian(X, F));
            Xold = std::move(X);
            X = Xold - JacInv * F;
        }

        if (k % 10 == 0)
        {
            JacInv = JacInvNew = math::pInv(problem.jacobian(X, F));
            X = Xold - JacInv * F;
        }

        if (blaze::linfNorm(F) <= problem.tolerance())
            found = true;
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
    double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
    std::size_t j = 0UL, k = 0UL;
    constexpr std::size_t k_max = 300UL;
    std::vector<double> h_k;
    h_k.reserve(k_max);

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
            if (++isfinite_guard > 100UL)
                break;
        }

        D_inv = math::pInv(D);
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
            h_k.push_back(h);

            if (h <= (h_0 - improvementFactor))
                break;
            else
                j++;
        }

        lambda = blaze::pow(0.5, static_cast<double>(h_k.size() - 1UL));
        initGuess -= lambda * d;
        h_k.clear();
        j = 0UL;

        if (blaze::linfNorm(f) <= problem.tolerance())
            return makeResult(true, k, f);
    }

    if (!found)
    {
        sanitizeBVPGuess(initGuess);
        initGuess *= 0.75;
        FKResult fallback = PowellDogLegSolver{}.solve(initGuess, problem);

        if (!fallback)
        {
            sanitizeBVPGuess(initGuess);
            initGuess *= 0.75;
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
