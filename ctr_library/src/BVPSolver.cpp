// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include "BVPSolver.hpp"
#include "CTR.hpp"
#include "mathOperations.hpp"

#include <cmath>
#include <vector>

namespace ctr {

// ─── Utility shared by all solvers ───────────────────────────────────────────

// Use the shared sanitizer defined in CTRTypes.hpp to avoid duplication with CTR.cpp.
using detail::sanitizeBVPGuess;

// ─── Powell Dog-Leg ──────────────────────────────────────────────────────────

bool PowellDogLegSolver::solve(bvp_type &initGuess, CTR &ctr)
{
    bool found;
    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    double alpha, beta, delta, eps1, eps2, rho, c;
    bvp_type g, f, f_new, x_new, h_sd, h_gn, h_dl;
    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM> J;

    sanitizeBVPGuess(initGuess);

    delta = 1.0;
    eps1 = eps2 = 1.0e-22;

    f = ctr.ODESolver(initGuess);
    J = ctr.jac_BVP(initGuess, f);
    g = blaze::trans(J) * f;

    found = ((blaze::linfNorm(f) <= ctr.accuracy()) || (blaze::linfNorm(g) <= eps1));

    while (!found && (k < k_max))
    {
        k++;

        alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
        h_sd  = -alpha * g;
        h_gn  = -mathOp::pInv(J) * f;

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
                    beta = (-c + std::sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) / blaze::sqrNorm(h_gn - h_sd);
                else
                    beta = (delta * delta - blaze::sqrNorm(h_sd)) / (c + std::sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
                h_dl = h_sd + beta * (h_gn - h_sd);
            }
        }

        if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
            found = true;
        else
        {
            x_new = initGuess + h_dl;
            f_new = ctr.ODESolver(x_new);
            rho   = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h_dl) * ((delta * h_dl) - g));

            if (rho > 0.0)
            {
                initGuess = std::move(x_new);
                f         = std::move(f_new);
                J         = ctr.jac_BVP(initGuess, f);
                g         = blaze::trans(J) * f;
                if ((blaze::linfNorm(f) <= ctr.accuracy()) || (blaze::linfNorm(g) <= eps1))
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

    return found;
}

// ─── Levenberg-Marquardt ─────────────────────────────────────────────────────

bool LevenbergMarquardtSolver::solve(bvp_type &initGuess, CTR &ctr)
{
    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    bvp_type h, g, f, f_new;
    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::columnMajor> J, A;
    blaze::IdentityMatrix<double> I(BVP_DIM);
    double rho, nu = 2.0, mu, tau = 1.0e-3, e1 = 1.0e-18, e2 = 1.0e-25;
    bool found;

    sanitizeBVPGuess(initGuess);

    f    = ctr.ODESolver(initGuess);
    J    = ctr.jac_BVP(initGuess, f);
    A    = blaze::trans(J) * J;
    g    = blaze::trans(J) * f;
    found = (blaze::linfNorm(g) <= e1);
    mu   = tau * blaze::max(blaze::diagonal(A));

    while (!found && (k < k_max))
    {
        k++;
        blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

        f_new = ctr.ODESolver(initGuess + h);
        rho   = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h) * ((mu * h) - g));

        if (rho > 0.0)
        {
            initGuess += h;
            J     = ctr.jac_BVP(initGuess, f_new);
            A     = blaze::trans(J) * J;
            f     = std::move(f_new);
            g     = blaze::trans(J) * f;
            found = (blaze::linfNorm(g) <= e1);
            mu    = mu * std::max(0.33333333, 1.0 - blaze::pow(2.0 * rho - 1.0, 3.0));
            nu    = 2.0;
        }
        else
        {
            mu = mu * nu;
            nu = 2.0 * nu;
        }

        if (blaze::linfNorm(f) <= ctr.accuracy())
            found = true;
    }

    return found;
}

// ─── Broyden ─────────────────────────────────────────────────────────────────

bool BroydenSolver::solve(bvp_type &initGuess, CTR &ctr)
{
    bool found;
    sanitizeBVPGuess(initGuess);

    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::columnMajor> JacInv, JacInvNew;
    bvp_type F, Fold, X, Xold, deltaX, deltaF;

    F         = ctr.ODESolver(initGuess);
    X         = std::move(initGuess);
    JacInvNew = JacInv = mathOp::pInv(ctr.jac_BVP(X, F));
    found     = (blaze::linfNorm(F) <= ctr.accuracy());

    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    while (!found && (k < k_max))
    {
        k++;

        deltaX = X - Xold;
        deltaF = F - Fold;

        JacInv = std::move(JacInvNew);
        if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.0))
            JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) * blaze::trans(deltaX) * JacInv;
        else
            JacInvNew = JacInv;

        Xold = std::move(X);
        Fold = std::move(F);
        X    = Xold - JacInv * F;
        F    = ctr.ODESolver(X);

        while (blaze::isnan(F))
        {
            X *= 0.75;
            sanitizeBVPGuess(X);
            F         = ctr.ODESolver(X);
            JacInv    = JacInvNew = mathOp::pInv(ctr.jac_BVP(X, F));
            Xold      = std::move(X);
            X         = Xold - JacInv * F;
        }

        if (k % 10 == 0)
        {
            JacInv = JacInvNew = mathOp::pInv(ctr.jac_BVP(X, F));
            X = Xold - JacInv * F;
        }

        if (blaze::linfNorm(F) <= ctr.accuracy())
            found = true;
    }

    initGuess = std::move(X);
    return found;
}

// ─── Broyden II ──────────────────────────────────────────────────────────────

bool BroydenIISolver::solve(bvp_type &initGuess, CTR &ctr)
{
    bool found;
    sanitizeBVPGuess(initGuess);

    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::columnMajor> Jac, JacNew;
    bvp_type F, Fold, X, Xold, deltaX, deltaF;

    F      = ctr.ODESolver(initGuess);
    X      = std::move(initGuess);
    JacNew = ctr.jac_BVP(initGuess, F);
    found  = (blaze::linfNorm(F) <= ctr.accuracy());

    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;
    while (!found && (k < k_max))
    {
        k++;

        deltaX = X - Xold;
        deltaF = F - Fold;

        Jac = std::move(JacNew);
        if (blaze::sqrNorm(deltaX) > 0.0)
            JacNew = Jac + blaze::sqrNorm(deltaX) * (deltaF - (Jac * deltaX)) * blaze::trans(deltaX);
        else
            JacNew = Jac;

        Xold = std::move(X);
        Fold = std::move(F);
        X    = Xold - mathOp::pInv(Jac) * F;
        F    = ctr.ODESolver(X);

        while (blaze::isnan(F))
        {
            X *= 0.75;
            sanitizeBVPGuess(X);
            F      = ctr.ODESolver(X);
            JacNew = ctr.jac_BVP(X, F);
            Xold   = std::move(X);
            X      = Xold - mathOp::pInv(Jac) * F;
        }

        if (k % 10 == 0)
        {
            JacNew = ctr.jac_BVP(X, F);
            X = Xold - mathOp::pInv(Jac) * F;
        }

        if (blaze::linfNorm(F) <= ctr.accuracy())
            found = true;
    }

    initGuess = std::move(X);
    return found;
}

// ─── Newton-Raphson ──────────────────────────────────────────────────────────

bool NewtonRaphsonSolver::solve(bvp_type &initGuess, CTR &ctr)
{
    bool found;
    bvp_type Residue, Residue_new, d_Residue, int_Residue, dGuess;

    sanitizeBVPGuess(initGuess);

    Residue = ctr.ODESolver(initGuess);
    found   = (blaze::linfNorm(Residue) <= ctr.accuracy());

    blaze::DiagonalMatrix<blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::rowMajor>> Kp, Ki, Kd;
    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::columnMajor> jac_bvp;
    blaze::diagonal(Kp) = 0.450;
    blaze::diagonal(Ki) = 0.005;
    blaze::diagonal(Kd) = 0.002;

    std::size_t k = 0UL;
    constexpr std::size_t k_max = 300UL;

    while (!found && (k < k_max))
    {
        k++;
        jac_bvp   = ctr.jac_BVP(initGuess, Residue);
        dGuess    = mathOp::pInv(jac_bvp) * (Kp * Residue + Ki * int_Residue + Kd * d_Residue);
        initGuess -= dGuess;

        sanitizeBVPGuess(initGuess);

        Residue_new  = ctr.ODESolver(initGuess);
        d_Residue    = Residue_new - Residue;
        int_Residue += Residue_new;
        Residue      = std::move(Residue_new);

        if (blaze::linfNorm(Residue) <= ctr.accuracy())
            found = true;
    }

    return found;
}

// ─── Modified Newton-Raphson ─────────────────────────────────────────────────

bool ModifiedNewtonRaphsonSolver::solve(bvp_type &initGuess, CTR &ctr)
{
    /*
     * Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition
     * by Josef Stoer & Roland Bulirsch.
     */
    sanitizeBVPGuess(initGuess);

    bool found;
    bvp_type f(ctr.ODESolver(initGuess)), d;
    blaze::StaticVector<double, BVP_DIM, blaze::rowVector> Dh;
    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM> D, D_inv;
    double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
    std::size_t j = 0UL, k = 0UL;
    constexpr std::size_t k_max = 300UL;
    std::vector<double> h_k;
    h_k.reserve(k_max);

    found = (blaze::linfNorm(f) <= ctr.accuracy());

    auto setupMethod = [&]() -> void
    {
        f = ctr.ODESolver(initGuess);
        D = ctr.jac_BVP(initGuess, f);

        std::size_t isfinite_guard = 0UL;
        while (!blaze::isfinite(D))
        {
            initGuess *= 0.75;
            sanitizeBVPGuess(initGuess);
            f = ctr.ODESolver(initGuess);
            D = ctr.jac_BVP(initGuess, f);
            if (++isfinite_guard > 100UL) break;
        }

        D_inv   = mathOp::pInv(D);
        d       = D_inv * f;
        gamma   = 1.0 / (blaze::norm(D_inv) * blaze::norm(D));
        h_0     = blaze::sqrNorm(f);
        Dh      = 2.0 * blaze::trans(f) * D;
        d_norm  = blaze::norm(d);
        Dh_norm = blaze::norm(Dh);
    };

    while (!found && (k < k_max))
    {
        k++;
        setupMethod();

        while (true)
        {
            f = ctr.ODESolver(initGuess - blaze::pow(0.5, j) * d);

            while (!blaze::isfinite(f))
            {
                j++;
                f = ctr.ODESolver(initGuess - blaze::pow(0.5, j) * d);
                if (j > 10UL)
                {
                    initGuess *= 0.75;
                    sanitizeBVPGuess(initGuess);
                    setupMethod();
                    j = 0UL;
                }
            }

            h                 = blaze::sqrNorm(f);
            improvementFactor = blaze::pow(0.5, j) * 0.25 * gamma * d_norm * Dh_norm;
            h_k.push_back(h);

            if (h <= (h_0 - improvementFactor))
                break;
            else
                j++;
        }

        lambda     = blaze::pow(0.5, static_cast<double>(h_k.size() - 1UL));
        initGuess -= lambda * d;
        h_k.clear();
        j = 0UL;

        if (blaze::linfNorm(f) <= ctr.accuracy())
            return true;
    }

    if (!found)
    {
        sanitizeBVPGuess(initGuess);
        initGuess *= 0.75;
        found = PowellDogLegSolver{}.solve(initGuess, ctr);

        if (!found)
        {
            sanitizeBVPGuess(initGuess);
            initGuess *= 0.75;
            found = LevenbergMarquardtSolver{}.solve(initGuess, ctr);
        }
    }

    return found;
}

// ─── Factory ─────────────────────────────────────────────────────────────────

std::unique_ptr<BVPSolver> makeBVPSolver(mathOp::rootFindingMethod method)
{
    switch (method)
    {
    case mathOp::rootFindingMethod::NEWTON_RAPHSON:
        return std::make_unique<NewtonRaphsonSolver>();
    case mathOp::rootFindingMethod::LEVENBERG_MARQUARDT:
        return std::make_unique<LevenbergMarquardtSolver>();
    case mathOp::rootFindingMethod::POWELL_DOG_LEG:
        return std::make_unique<PowellDogLegSolver>();
    case mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON:
        return std::make_unique<ModifiedNewtonRaphsonSolver>();
    case mathOp::rootFindingMethod::BROYDEN:
        return std::make_unique<BroydenSolver>();
    case mathOp::rootFindingMethod::BROYDEN_II:
        return std::make_unique<BroydenIISolver>();
    }
    return std::make_unique<PowellDogLegSolver>(); // unreachable, appeases the compiler
}

} // namespace ctr
