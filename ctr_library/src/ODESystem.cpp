// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "ODESystem.hpp"

// default constructor
ODESystem::ODESystem() : m_u_ast_x(0.00), m_u_ast_y(0.00), m_EI(0.00), m_GJ(0.00), m_f(0.00)
{
	m_e3 = {0.00, 0.00, 1.00};
}

// overloaded constructor
ODESystem::ODESystem(const blaze::StaticVector<double, 3UL> &u_ast_x, const blaze::StaticVector<double, 3UL> &u_ast_y, const blaze::StaticVector<double, 3UL> &EI, const blaze::StaticVector<double, 3UL> &GJ) : m_u_ast_x(u_ast_x), m_u_ast_y(u_ast_y), m_EI(EI), m_GJ(GJ)
{
	m_e3 = {0.00, 0.00, 1.00};
	m_f = 0.00;
}

// copy constructor
ODESystem::ODESystem(const ODESystem &rhs) : m_u_ast_x(rhs.m_u_ast_x), m_u_ast_y(rhs.m_u_ast_y), m_EI(rhs.m_EI), m_GJ(rhs.m_GJ), m_e3(rhs.m_e3), m_f(rhs.m_f) {}

// move constructor
ODESystem::ODESystem(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_f = std::move(rhs.m_f);
	}
}

// ODESystem destructor
ODESystem::~ODESystem()
{
	// nothing to be done
}

// Copy assignment operator
ODESystem &ODESystem::operator=(const ODESystem &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = rhs.m_u_ast_x;
		this->m_u_ast_y = rhs.m_u_ast_y;
		this->m_EI = rhs.m_EI;
		this->m_GJ = rhs.m_GJ;
		this->m_e3 = rhs.m_e3;
		this->m_f = rhs.m_f;
	}
	return *this;
}

// move assignment operator
ODESystem &ODESystem::operator=(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_f = std::move(rhs.m_f);
	}
	return *this;
}

// functor that implements the system of ODEs governing a three-tube CTR
void ODESystem::operator()(const state_type &y, state_type &dyds, const double s)
{
	// 1st element of y computes the bending moment of the first (innermost) tube along the x direction
	// 2nd element of y computes the bending moment of the first (innermost) tube along the y direction
	// next 3 elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
	// next 2 elements of y are twist angles, theta_i = [theta_1 theta_2  theta_3]
	// last 7 elements are r(position) and h(quaternion-orientations) of the local frame, respectively at each arc-length s

	double dtheta_2 = y[3UL] - y[2UL];
	double dtheta_3 = y[4UL] - y[2UL];

	// implementing curvature equation u_i = transpose(R_z(theta_i))*u_1 + \dot{theta_i}*e3
	blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1, R2(mathOp::rotz(y[6UL])), R3(mathOp::rotz(y[7UL]));
	blaze::StaticVector<double, 3UL> u1, u2, u3, mb = {y[0UL], y[1UL], 0.00};

	blaze::StaticVector<double, 3UL> u1_ast = {m_u_ast_x[0UL], m_u_ast_y[0UL], 0.00};
	blaze::StaticVector<double, 3UL> u2_ast = {m_u_ast_x[1UL], m_u_ast_y[1UL], 0.00};
	blaze::StaticVector<double, 3UL> u3_ast = {m_u_ast_x[2UL], m_u_ast_y[2UL], 0.00};

	// estimating curvature of the first tube along the x and y directions
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K_inv, K1, K2, K3;

	K1(0UL, 0UL) = K1(1UL, 1UL) = m_EI[0UL];
	K1(2UL, 2UL) = m_GJ[0UL];
	K2(0UL, 0UL) = K2(1UL, 1UL) = m_EI[1UL];
	K2(2UL, 2UL) = m_GJ[1UL];
	K3(0UL, 0UL) = K3(1UL, 1UL) = m_EI[2UL];
	K3(2UL, 2UL) = m_GJ[2UL];

	K_inv = blaze::inv(K1 + K2 + K3);

	mb[2UL] = K1(2UL, 2UL) * y[2UL] + K2(2UL, 2UL) * y[3UL] + K3(2UL, 2UL) * y[4UL];
	// estimating the curvature of the innermost tube along the x,y directions
	u1 = K_inv * (mb + (K1 * u1_ast) + (R2 * K2 * u2_ast) + (R3 * K3 * u3_ast)); // expression valid only along the x,y directions
	// grabbing the torsion along the z-direction from state vector
	u1[2UL] = y[2UL];

	// curvatures of the intermediate and outermost tubes
	u2 = blaze::trans(R2) * u1 + (dtheta_2 * m_e3);
	u3 = blaze::trans(R3) * u1 + (dtheta_3 * m_e3);

	// gets orientation of the innermost tube (Tb 1) at the current arc-length
	mathOp::getSO3(blaze::subvector<11UL, 4UL>(y), R1);

	// estimating the twist curvatures (uz_i) and twist angles (theta_i)
	auto computeTwists = [&](size_t idx, const blaze::StaticVector<double, 3UL> &u) -> void
	{
		if (m_GJ[idx] != 0.00)
		{
			// uz_i = ( (E_i * I_i) / (G_i * J_i) ) * (ux_i * uy_ast - uy_i * ux_ast)
			dyds[2 + idx] = (m_EI[idx] / m_GJ[idx]) * (u[0UL] * m_u_ast_y[idx] - u[1UL] * m_u_ast_x[idx]);
			// dtheta_i = uz_i - uz_1
			dyds[5 + idx] = u[2UL] - u1[2UL];
		}
		else
			dyds[2 + idx] = dyds[5 + idx] = 0.00;
	};

	computeTwists(0UL, u1);
	computeTwists(1UL, u2);
	computeTwists(2UL, u3);

	// internal moment of tube 1 along the x and y directions
	blaze::subvector<0UL, 2UL>(dyds) = blaze::subvector<0UL, 2UL>(-mathOp::hatOperator(u1) * mb - mathOp::hatPreMultiply(m_e3, blaze::trans(R1)) * this->m_f); // - u^mb - e3^R1'F

	// spatial derivative of the quaternion representation h_dot
	blaze::subvector<11UL, 4UL>(dyds) = mathOp::quaternionDiff(u1, blaze::subvector<11UL, 4UL>(y));

	// calculating r_dot
	blaze::subvector<8UL, 3UL>(dyds) = blaze::column<2UL>(R1); // r_dot = R1 * e3
}

// setter method for updating the parameters forward kinematics computation
void ODESystem::setEquationParameters(const blaze::StaticVector<double, 3UL> &u_ast_x,
									  const blaze::StaticVector<double, 3UL> &u_ast_y,
									  const blaze::StaticVector<double, 3UL> &EI,
									  const blaze::StaticVector<double, 3UL> &GJ,
									  const blaze::StaticVector<double, 3UL> &force)
{
	this->m_u_ast_x = u_ast_x;
	this->m_u_ast_y = u_ast_y;
	this->m_EI = EI;
	this->m_GJ = GJ;
	this->m_f = force;
}