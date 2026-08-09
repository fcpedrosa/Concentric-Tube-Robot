#pragma once

#include <blaze/Math.h>
#include "ctr/Types.hpp"
#include "ctr/Segment.hpp"

namespace ctr
{

/**
 * @brief Implements the system of Ordinary Differential Equations (ODEs) that model the
 *        kinematics of a three-tube CTR.
 *
 * Designed as a callable functor compatible with Boost.Odeint integrators.
 * Reparameterized via setEquationParameters() before integrating each
 * arc-length segment; all quantities that are constant within a segment
 * (stiffness sums, EI·u* products, EI/GJ ratios) are precomputed there, so
 * operator() — evaluated four times per integration step — is pure scalar
 * arithmetic with two sin/cos pairs and no matrix algebra.
 *
 * @note Advanced API: most users never touch this class directly — CTR owns
 *       and drives it internally.
 */
class ODESystem
{
  private:
    // ─── Raw per-segment inputs ───────────────────────────────────────────────
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_x; ///< Pre-curvature x per tube [1/m]; 0 if absent/straight.
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_y; ///< Pre-curvature y per tube [1/m].
    blaze::StaticVector<double, NUM_TUBES> m_GJ;      ///< Torsional stiffness per tube [N·m²]; 0 if absent.

    // ─── Per-segment invariants (precomputed in setEquationParameters) ────────
    blaze::StaticVector<double, NUM_TUBES> m_EIoverGJ; ///< EIᵢ/GJᵢ, 0 where the tube is absent.
    blaze::StaticVector<double, NUM_TUBES> m_EIux;     ///< EIᵢ·u*ₓᵢ [N·m].
    blaze::StaticVector<double, NUM_TUBES> m_EIuy;     ///< EIᵢ·u*ᵧᵢ [N·m].
    double m_invSumEI{0.0};                            ///< 1 / ΣEIᵢ [1/(N·m²)].

    blaze::StaticVector<double, 3UL> m_f; ///< Distal point force [N], global frame.
    bool m_hasForce{false};               ///< Short-circuits the force term when unloaded.

  public:
    /** @brief Unit vector in the z-direction (constant). Exposed for use in boundary conditions. */
    static constexpr blaze::StaticVector<double, 3UL> kE3{0.0, 0.0, 1.0};

    /** @brief Default constructor. Initialises all parameters to zero. */
    ODESystem();

    ODESystem(const ODESystem &) = default;     ///< Copyable.
    ODESystem(ODESystem &&) noexcept = default; ///< Movable.
    ~ODESystem() = default;
    ODESystem &operator=(const ODESystem &) = default;     ///< Copy-assignable.
    ODESystem &operator=(ODESystem &&) noexcept = default; ///< Move-assignable.

    /**
     * @brief Functor implementing the ODE right-hand side for Boost.Odeint.
     *
     * @param y    Current state vector at arc-length s.
     * @param dyds Output: spatial derivative of the state vector at arc-length s.
     * @param s    Arc-length (not used explicitly but required by Boost.Odeint interface).
     */
    void operator()(const state_type &y, state_type &dyds, double s) const noexcept;

    /**
     * @brief Loads one segment's parameters and precomputes all per-segment invariants.
     *
     * @param seg    The segment decomposition of the backbone.
     * @param segIdx Index of the segment about to be integrated.
     * @param force  Distal-end point force in the global frame [N].
     */
    void setEquationParameters(const Segment &seg, std::size_t segIdx, const blaze::StaticVector<double, 3UL> &force);
};

} // namespace ctr
