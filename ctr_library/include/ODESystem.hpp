#pragma once

#include <blaze/Math.h>
#include "CTRTypes.hpp"
#include "mathOperations.hpp"

namespace ctr
{

/**
 * @brief Implements the system of Ordinary Differential Equations (ODEs) that model the
 *        kinematics of a three-tube CTR.
 *
 * Designed as a callable functor compatible with Boost.Odeint integrators.
 */
class ODESystem
{
  private:
    /**
     * @brief Pre-curvature of the tubes along the 'x' direction in the current segment.
     */
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_x;

    /**
     * @brief Pre-curvature of the tubes along the 'y' direction in the current segment.
     */
    blaze::StaticVector<double, NUM_TUBES> m_u_ast_y;

    /**
     * @brief Bending stiffness of each tube in the current segment. Zero if the tube is absent.
     */
    blaze::StaticVector<double, NUM_TUBES> m_EI;

    /**
     * @brief Torsional stiffness of each tube in the current segment. Zero if the tube is absent.
     */
    blaze::StaticVector<double, NUM_TUBES> m_GJ;

    /**
     * @brief Point force acting at the distal end of the CTR.
     */
    blaze::StaticVector<double, 3UL> m_f;

  public:
    /** @brief Unit vector in the z-direction (constant). Exposed for use in boundary conditions. */
    static constexpr blaze::StaticVector<double, 3UL> kE3{0.0, 0.0, 1.0};
    /**
     * @brief Default constructor. Initialises all parameters to zero.
     */
    ODESystem();

    /**
     * @brief Overloaded constructor.
     *
     * @param u_ast_x Pre-curvature along x for each tube in the current segment.
     * @param u_ast_y Pre-curvature along y for each tube in the current segment.
     * @param EI      Bending stiffness for each tube in the current segment.
     * @param GJ      Torsional stiffness for each tube in the current segment.
     */
    ODESystem(const blaze::StaticVector<double, NUM_TUBES> &u_ast_x,
              const blaze::StaticVector<double, NUM_TUBES> &u_ast_y, const blaze::StaticVector<double, NUM_TUBES> &EI,
              const blaze::StaticVector<double, NUM_TUBES> &GJ);

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
    void operator()(const state_type &y, state_type &dyds, double s) noexcept;

    /**
     * @brief Updates all kinematic parameters before integrating a new segment.
     *
     * @param u_ast_x Pre-curvature along x for each tube in the new segment.
     * @param u_ast_y Pre-curvature along y for each tube in the new segment.
     * @param EI      Bending stiffness for each tube in the new segment.
     * @param GJ      Torsional stiffness for each tube in the new segment.
     * @param force   Distal-end point force in the global frame.
     */
    void setEquationParameters(const blaze::StaticVector<double, NUM_TUBES> &u_ast_x,
                               const blaze::StaticVector<double, NUM_TUBES> &u_ast_y,
                               const blaze::StaticVector<double, NUM_TUBES> &EI,
                               const blaze::StaticVector<double, NUM_TUBES> &GJ,
                               const blaze::StaticVector<double, 3UL> &force);
};

} // namespace ctr
