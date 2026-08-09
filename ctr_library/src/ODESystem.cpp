// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "ODESystem.hpp"

namespace ctr {

// kE3 out-of-line definition required for ODR in C++17; C++20 constexpr inline makes it implicit
// constexpr blaze::StaticVector<double, 3UL> ODESystem::kE3;  // uncomment if linking errors on older ABIs

// default constructor
ODESystem::ODESystem()
	: m_u_ast_x(0.0), m_u_ast_y(0.0), m_EI(0.0), m_GJ(0.0), m_f(0.0)
{}

// overloaded constructor
ODESystem::ODESystem(const blaze::StaticVector<double, NUM_TUBES> &u_ast_x,
					 const blaze::StaticVector<double, NUM_TUBES> &u_ast_y,
					 const blaze::StaticVector<double, NUM_TUBES> &EI,
					 const blaze::StaticVector<double, NUM_TUBES> &GJ)
	: m_u_ast_x(u_ast_x), m_u_ast_y(u_ast_y), m_EI(EI), m_GJ(GJ), m_f(0.0)
{}

// functor that implements the system of ODEs governing a three-tube CTR
void ODESystem::operator()(const state_type &y, state_type &dyds, const double /*s*/) noexcept
{
	using namespace StateIdx;

	// 1st element of y computes the bending moment of the first (innermost) tube along the x direction
	// 2nd element of y computes the bending moment of the first (innermost) tube along the y direction
	// next 3 elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
	// next 2 elements of y are twist angles, theta_i = [theta_1 theta_2  theta_3]
	// last 7 elements are r(position) and h(quaternion-orientations) of the local frame, respectively at each arc-length s

	// dθᵢ/ds = u_iz − u_1z  (Rucker et al. 2010, eq. 16): torsional-curvature difference,
	// NOT a difference of the twist angles themselves.
	const double dtheta_2 = y[UZ_2] - y[UZ_1];
	const double dtheta_3 = y[UZ_3] - y[UZ_1];

	// implementing curvature equation u_i = transpose(R_z(theta_i))*u_1 + \dot{theta_i}*e3
	blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1;
	const blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R2(mathOp::rotz(y[THETA_2]));
	const blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R3(mathOp::rotz(y[THETA_3]));

	const blaze::StaticVector<double, 3UL> u1_ast = {m_u_ast_x[0UL], m_u_ast_y[0UL], 0.0};
	const blaze::StaticVector<double, 3UL> u2_ast = {m_u_ast_x[1UL], m_u_ast_y[1UL], 0.0};
	const blaze::StaticVector<double, 3UL> u3_ast = {m_u_ast_x[2UL], m_u_ast_y[2UL], 0.0};

	// estimating curvature of the first tube along the x and y directions
	const blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K1 = {
		{m_EI[0UL], 0.0, 0.0},
		{0.0, m_EI[0UL], 0.0},
		{0.0, 0.0, m_GJ[0UL]}};

	const blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K2 = {
		{m_EI[1UL], 0.0, 0.0},
		{0.0, m_EI[1UL], 0.0},
		{0.0, 0.0, m_GJ[1UL]}};

	const blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K3 = {
		{m_EI[2UL], 0.0, 0.0},
		{0.0, m_EI[2UL], 0.0},
		{0.0, 0.0, m_GJ[2UL]}};

	const blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K_inv{blaze::inv(K1 + K2 + K3)};

	const blaze::StaticVector<double, 3UL> mb{y[MB_X],
											  y[MB_Y],
											  K1(2UL, 2UL) * y[UZ_1] + K2(2UL, 2UL) * y[UZ_2] + K3(2UL, 2UL) * y[UZ_3]};

	// estimating the curvature of the innermost tube along the x,y directions
	blaze::StaticVector<double, 3UL> u1 = K_inv * (mb + (K1 * u1_ast) + (R2 * K2 * u2_ast) + (R3 * K3 * u3_ast));
	// grabbing the torsion along the z-direction from state vector
	u1[2UL] = y[UZ_1];

	// curvatures of the intermediate and outermost tubes
	const blaze::StaticVector<double, 3UL> u2 = blaze::trans(R2) * u1 + (dtheta_2 * kE3);
	const blaze::StaticVector<double, 3UL> u3 = blaze::trans(R3) * u1 + (dtheta_3 * kE3);

	// gets orientation of the innermost tube (Tb 1) at the current arc-length
	mathOp::getSO3(blaze::subvector<QUAT_W, 4UL>(y), R1);

	// estimating the twist curvatures (uz_i) and twist angles (theta_i)
	auto computeTwists = [&](std::size_t idx, const blaze::StaticVector<double, 3UL> &u) -> void
	{
		if (m_GJ[idx] != 0.0)
		{
			// uz_i = ( (E_i * I_i) / (G_i * J_i) ) * (ux_i * uy_ast - uy_i * ux_ast)
			dyds[UZ_1 + idx]    = (m_EI[idx] / m_GJ[idx]) * (u[0UL] * m_u_ast_y[idx] - u[1UL] * m_u_ast_x[idx]);
			// dtheta_i = uz_i - uz_1
			dyds[THETA_1 + idx] = u[2UL] - u1[2UL];
		}
		else
		{
			dyds[UZ_1 + idx] = dyds[THETA_1 + idx] = 0.0;
		}
	};

	computeTwists(0UL, u1);
	computeTwists(1UL, u2);
	computeTwists(2UL, u3);

	// internal moment of tube 1 along the x and y directions
	blaze::subvector<MB_X, 2UL>(dyds) =
		blaze::subvector<0UL, 2UL>(-mathOp::hatOperator(u1) * mb - mathOp::hatPreMultiply(kE3, blaze::trans(R1)) * m_f);

	// spatial derivative of the quaternion representation h_dot
	blaze::subvector<QUAT_W, 4UL>(dyds) = mathOp::quaternionDiff(u1, blaze::subvector<QUAT_W, 4UL>(y));

	// calculating r_dot = R1 * e3
	blaze::subvector<POS_X, 3UL>(dyds) = blaze::column<2UL>(R1);
}

void ODESystem::setEquationParameters(const blaze::StaticVector<double, NUM_TUBES> &u_ast_x,
									  const blaze::StaticVector<double, NUM_TUBES> &u_ast_y,
									  const blaze::StaticVector<double, NUM_TUBES> &EI,
									  const blaze::StaticVector<double, NUM_TUBES> &GJ,
									  const blaze::StaticVector<double, 3UL>       &force)
{
	m_u_ast_x = u_ast_x;
	m_u_ast_y = u_ast_y;
	m_EI      = EI;
	m_GJ      = GJ;
	m_f       = force;
}

} // namespace ctr
