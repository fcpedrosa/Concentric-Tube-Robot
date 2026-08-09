#include "ctr/detail/mathOperations.hpp"
#include <iostream>

namespace ctr::math
{

/**
 * @brief Damped least-squares pseudo-inverse via SVD.
 *
 * Computes the Moore-Penrose pseudo-inverse of M using the decomposition
 *   M = U × diag(s) × Vᵀ
 * and the damped inversion formula
 *   σ_inv = σ / (σ² + λ²),  λ² = 1×10⁻²⁵.
 *
 * Placing the full SVD body here (rather than in the header) prevents the
 * expensive Blaze SVD instantiation from being compiled into every translation
 * unit that includes mathOperations.hpp.
 */
blaze::DynamicMatrix<double, blaze::columnMajor> pInv(blaze::DynamicMatrix<double, blaze::columnMajor> Mcm)
{
    blaze::DynamicMatrix<double, blaze::columnMajor> U, V;
    blaze::DynamicVector<double> s;

    try
    {
        blaze::svd(Mcm, U, s, V);
    }
    catch (const std::exception &e)
    {
        std::cerr << "ctr::math::pInv — Blaze SVD failed: " << e.what() << '\n';
        // Return a zero matrix of the transposed shape so callers do not crash.
        return blaze::DynamicMatrix<double, blaze::columnMajor>(Mcm.columns(), Mcm.rows(), 0.0);
    }

    constexpr double lambda_sq = 1.0E-25; // damping term
    const std::size_t n = s.size();
    blaze::DynamicMatrix<double, blaze::columnMajor> S_inv(U.columns(), V.rows(), 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double d = s[i];
        S_inv(i, i) = (d <= 1.0E-25) ? 0.0 : d / (d * d + lambda_sq);
    }

    return blaze::trans(blaze::evaluate(U * S_inv * V));
}

} // namespace ctr::math
