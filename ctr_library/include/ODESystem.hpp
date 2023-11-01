#pragma once

#include <blaze/Math.h>
#include "mathOperations.hpp"

typedef blaze::StaticVector<double, 15UL> state_type;

class ODESystem
{
private:
	blaze::StaticVector<double, 3UL> m_u_ast_x, m_u_ast_y, m_EI, m_GJ, m_e3, m_f;
	// matrices for the derivative propagation approach

public:
	// default constructor
	ODESystem();

	// overloaded constructor
	ODESystem(const blaze::StaticVector<double, 3UL> &u_ast_x, const blaze::StaticVector<double, 3UL> &u_ast_y, const blaze::StaticVector<double, 3UL> &EI, const blaze::StaticVector<double, 3UL> &GJ);

	// copy constructor
	ODESystem(const ODESystem &rhs);

	// move constructor
	ODESystem(ODESystem &&rhs) noexcept;

	// ODESystem destructor
	~ODESystem() = default;

	// Copy assignment operator
	ODESystem &operator=(const ODESystem &rhs);

	// move assignment operator
	ODESystem &operator=(ODESystem &&rhs) noexcept;

	// functor that implements the system of ODEs governing a three-tube CTR
	void operator()(const state_type &y, state_type &dyds, const double s);

	// setter method for updating the parameters for forward kinematics computation
	void setEquationParameters(const blaze::StaticVector<double, 3UL> &u_ast_x,
							   const blaze::StaticVector<double, 3UL> &u_ast_y,
							   const blaze::StaticVector<double, 3UL> &EI,
							   const blaze::StaticVector<double, 3UL> &GJ,
							   const blaze::StaticVector<double, 3UL> &force);
};