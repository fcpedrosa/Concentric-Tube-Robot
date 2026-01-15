#pragma once

#include <blaze/Math.h>
#include "mathOperations.hpp"

typedef blaze::StaticVector<double, 15UL> state_type;

/**
 * @file ODESystem.hpp
 * @brief ODE system definitions for CTR kinematics.
 * @ingroup ctr_core
 * @details Provides the state dynamics used by the numerical integrator.
 */

/**
 * @brief Implements system of Ordinary Differential equations (ODEs) that model the kinematics of a three-tube CTR
 * @ingroup ctr_core
 */
class ODESystem
{
private:
	/**
	 * @brief Pre-curvature of the tubes along the 'x' direction in the current segment.
	 * The i-th entry of the vector corresponds to the pre-curvature in the x-direction for the i-th tube in the concentric arrangement.
	 */
	blaze::StaticVector<double, 3UL> m_u_ast_x;

	/**
	 * @brief Pre-curvature of the tubes along the 'y' direction in the current segment.
	 * The i-th entry of the vector corresponds to the pre-curvature in the y-direction for the i-th tube in the concentric arrangement.
	 */
	blaze::StaticVector<double, 3UL> m_u_ast_y;

	/**
	 * @brief Bending stiffness of each of the tubes in the current segment.
	 * If a tube is not present in the current segment, the corresponding entry will be zero.
	 */
	blaze::StaticVector<double, 3UL> m_EI;

	/**
	 * @brief Torsional stiffness of each tube in the current segment.
	 * If a tube is not present in the current segment, the corresponding entry will be zero.
	 */
	blaze::StaticVector<double, 3UL> m_GJ;

	/**
	 * @brief A unit vector in the z-direction.
	 */
	blaze::StaticVector<double, 3UL> m_e3;

	/**
	 * @brief Point force acting at the distal end of the CTR.
	 */
	blaze::StaticVector<double, 3UL> m_f;

public:
	/**
	 * @brief Implements the default constructor for the ODESystem class
	 */
	ODESystem();

	/**
	 * @brief Implements the overloaded constructor for the ODESystem class.
	 *
	 * @param u_ast_x A 3-dimensional static Blaze vector of the pre-curvature of the tubes along the 'x' direction in the present segment.
	 * @param u_ast_y A 3-dimensional static Blaze vector of the pre-curvature of the tubes along the 'x' direction in the present segment.
	 * @param EI A 3-dimensional static Blaze vector of bending stiffness of each one of the tubes in the present segment. If the i-th tube isn't present in the current segment, the i-th entry of EI_i will be zero.
	 * @param GJ A 3-dimensional static Blaze vector of torsional stiffness of each one of the tubes in the present segment. If the i-th tube isn't present in the current segment, the i-th entry of GJ_i will be zero.
	 */
	ODESystem(const blaze::StaticVector<double, 3UL> &u_ast_x, const blaze::StaticVector<double, 3UL> &u_ast_y, const blaze::StaticVector<double, 3UL> &EI, const blaze::StaticVector<double, 3UL> &GJ);

	/**
	 * @brief Implements the copy constructor for the ODESystem class.
	 *
	 * @param rhs The source ODESystem object to copy from.
	 */
	ODESystem(const ODESystem &rhs);

	/**
	 * @brief Implements the move constructor for the ODESystem class.
	 *
	 * @param rhs The source ODESystem object to move from.
	 */
	ODESystem(ODESystem &&rhs) noexcept;

	/**
	 * @brief Destroys the ODESystem object.
	 */
	~ODESystem() = default;

	/**
	 * @brief Implements the copy assignment operator for the ODESystem class.
	 *
	 * @param rhs The source ODESystem object to copy from.
	 * @return A reference to the assigned ODESystem object.
	 */
	ODESystem &operator=(const ODESystem &rhs);

	/**
	 * @brief Implements the move assignment operator for the ODESystem class.
	 *
	 * @param rhs The source ODESystem object to move from.
	 * @return A reference to the assigned ODESystem object.
	 */
	ODESystem &operator=(ODESystem &&rhs) noexcept;

	/**
	 * @brief Functor that overloads the constructor's signature and implements the system of ODEs governing a three-tube CTR
	 *
	 * @param[in] y A 15-dimensional static Blaze vector containing the current values of the state vector at the arc-length 's' within the current segment.
	 * @param[out] dyds A 15-dimensional static Blaze vector to be computed by the functor. Once the functor is executed, 'dyds' corresponds to the spatial derivative of the state vector at the arc-length 's'.
	 * @param[in] s The nonnegative scalar corresponding to the arc-length along the CTR backbone at which the computations are taking place
	 */
	void operator()(const state_type &y, state_type &dyds, const double s) noexcept;

	/**
	 * @brief Implements a setter method for updating the kinematic parameters before computation of the ODEs at each arc-length 's'
	 *
	 * @param[in] u_ast_x A 3-dimensional static Blaze vector of the pre-curvature of the tubes along the 'x' direction in the present segment.
	 * @param[in] u_ast_y A 3-dimensional static Blaze vector of the pre-curvature of the tubes along the 'x' direction in the present segment.
	 * @param[in] EI A 3-dimensional static Blaze vector of bending stiffness of each one of the tubes in the present segment. If the i-th tube isn't present in the current segment, the i-th entry of EI_i will be zero.
	 * @param[in] GJ A 3-dimensional static Blaze vector of torsional stiffness of each one of the tubes in the present segment. If the i-th tube isn't present in the current segment, the i-th entry of GJ_i will be zero.
	 * @param[in] force A 3-dimensional static Blaze vector containing the point force acting at the distal-end of the CTR.
	 */
	void setEquationParameters(const blaze::StaticVector<double, 3UL> &u_ast_x,
							   const blaze::StaticVector<double, 3UL> &u_ast_y,
							   const blaze::StaticVector<double, 3UL> &EI,
							   const blaze::StaticVector<double, 3UL> &GJ,
							   const blaze::StaticVector<double, 3UL> &force);
};