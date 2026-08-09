#pragma once

#include <blaze/Math.h>
#include "ctr/Types.hpp"

namespace ctr
{

/**
 * @brief Implements the system of Ordinary Differential Equations (ODEs) that model the
 *        kinematics of a three-tube CTR.
 *
 * Designed as a callable functor compatible with Boost.Odeint integrators.
 * Reparameterized via setEquationParameters() before integrating each
 * arc-length segment.
 *
 * @note Advanced API: most users never touch this class directly — CTR owns
 *       and drives it internally.
 */
class ODESystem
{
  private:
    /** @brief Pre-curvature of the tubes along the 'x' direction in the current segment [1/m]. */
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_x;

    /** @brief Pre-curvature of the tubes along the 'y' direction in the current segment [1/m]. */
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_y;

    /** @brief Bending stiffness of each tube in the current segment [N·m²]. Zero if the tube is absent. */
    blaze::StaticVector<double, NUM_TUBES> m_EI;

    /** @brief Torsional stiffness of each tube in the current segment [N·m²]. Zero if the tube is absent. */
    blaze::StaticVector<double, NUM_TUBES> m_GJ;

    /** @brief Point force acting at the distal end of the CTR [N], global frame. */
    blaze::StaticVector<double, 3UL> m_f;

  public:
    /** @brief Unit vector in the z-direction (constant). Exposed for use in boundary conditions. */
    static constexpr blaze::StaticVector<double, 3UL> kE3{0.0, 0.0, 1.0};

    /** @brief Default constructor. Initialises all parameters to zero. */
    ODESystem();

    ODESystem(const ODESystem &) = default;
    ODESystem(ODESystem &&) noexcept = default;
    ~ODESystem() = default;
    ODESystem &operator=(const ODESystem &) = default;
    ODESystem &operator=(ODESystem &&) noexcept = default;

    /**
     * @brief Functor implementing the ODE right-hand side for Boost.Odeint.
     *
     * @param y    Current state vector at arc-length s.
     * @param dyds Output: spatial derivative of the state vector at arc-length s.
     * @param s    Arc-length (not used explicitly but required by Boost.Odeint interface).
     */
    void operator()(const state_type &y, state_type &dyds, double s) const noexcept;

    /**
     * @brief Updates all kinematic parameters before integrating a new segment.
     *
     * @param u_ast_x Pre-curvature along x for each tube in the new segment [1/m].
     * @param u_ast_y Pre-curvature along y for each tube in the new segment [1/m].
     * @param EI      Bending stiffness for each tube in the new segment [N·m²].
     * @param GJ      Torsional stiffness for each tube in the new segment [N·m²].
     * @param force   Distal-end point force in the global frame [N].
     */
    void setEquationParameters(const blaze::StaticVector<double, NUM_TUBES> &u_ast_x,
                               const blaze::StaticVector<double, NUM_TUBES> &u_ast_y,
                               const blaze::StaticVector<double, NUM_TUBES> &EI,
                               const blaze::StaticVector<double, NUM_TUBES> &GJ,
                               const blaze::StaticVector<double, 3UL> &force);
};

} // namespace ctr
