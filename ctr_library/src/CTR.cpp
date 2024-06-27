// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include "CTR.hpp"

// overloaded class constructor
CTR::CTR(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, blaze::StaticVector<double, 6UL> &q, double Tol, mathOp::rootFindingMethod method) : m_accuracy(Tol), m_method(method), m_Tubes(Tb), m_beta(blaze::subvector<0UL, 3UL>(q)), m_q(q)
{
	this->m_theta_0 = {0.00, q[4UL] - q[3UL], q[5UL] - q[4UL]};
	this->m_e3 = {0.00, 0.00, 1.00};
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_h_0 = {1.00, 0.00, 0.00, 0.00};
	this->m_wf = 0.00;
	this->m_wm = 0.00;

	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// copy constructor
CTR::CTR(const CTR &rhs) : m_accuracy(rhs.m_accuracy), m_method(rhs.m_method), m_Tubes(rhs.m_Tubes), m_beta(rhs.m_beta),
						   m_q(rhs.m_q), m_theta_0(rhs.m_theta_0), m_e3(rhs.m_e3), m_h_0(rhs.m_h_0),
						   m_y(rhs.m_y), m_s(rhs.m_s), m_wf(rhs.m_wf), m_wm(rhs.m_wm)
{
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// move constructor
CTR::CTR(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_wf = std::move(rhs.m_wf);
		this->m_wm = std::move(rhs.m_wm);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}
}

// copy assignment operator
CTR &CTR::operator=(const CTR &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = rhs.m_Tubes;
		this->m_beta = rhs.m_beta;
		this->m_q = rhs.m_q;
		this->m_theta_0 = rhs.m_theta_0;
		this->m_wf = rhs.m_wf;
		this->m_wm = rhs.m_wm;
		this->m_e3 = rhs.m_e3;
		this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
		this->m_h_0 = rhs.m_h_0;
		this->m_y = rhs.m_y;
		this->m_s = rhs.m_s;
		this->m_stateEquations = std::make_unique<ODESystem>();
		this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
	}

	return *this;
}

// move assignment operator
CTR &CTR::operator=(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_wf = std::move(rhs.m_wf);
		this->m_wm = std::move(rhs.m_wm);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}

	return *this;
}

// function that resets the initial parameters for the ODESolver
void CTR::reset(const blaze::StaticVector<double, 5UL> &initGuess)
{
	blaze::StaticVector<double, 3UL> uz_0 = {initGuess[2UL], initGuess[3UL], initGuess[4UL]};
	// alpha1_0 =  alpha_1 - beta_1 * uz_1(0)
	double alpha1_0 = this->m_q[3UL] - this->m_beta[0UL] * uz_0[0UL];

	// clearing the observer's containers
	if (!this->m_s.empty())
	{
		this->m_y.clear();
		this->m_y.reserve(1000UL);
		this->m_s.clear();
		this->m_s.reserve(1000UL);
	}

	// theta_i(0) = alpha_1 - alpha_i - (beta_i * uz_i(0) - beta_1 * uz_1(0))
	this->m_theta_0 = {0.00, this->m_q[4UL] - this->m_beta[1UL] * uz_0[1UL] - alpha1_0, this->m_q[5UL] - this->m_beta[2UL] * uz_0[2UL] - alpha1_0};

	// transforming proximal orientation to quaternion representation
	mathOp::euler2Quaternion(0.00, alpha1_0, 0.00, this->m_h_0);
}

// function that solves (integrates) the CTR ode (state) equations
blaze::StaticVector<double, 5UL> CTR::ODESolver(const blaze::StaticVector<double, 5UL> &initGuess)
{
	// initGuess = [mb_x(0) mb_y(0) u1_z(0) u2_z(0) u3_z(0)]  --> vector of initial guesses for solving the BVP
	this->reset(initGuess); // resets CTR parameters and variables for a new iteration of the ode-solver

	// retrieving the bending & torsional stiffness and precurvatures in all segments of the CTR in the current
	auto [EI, GJ, U_x, U_y, S] = this->m_segment->returnParameters();

	// ##################################################### NUMERICAL METHODS FOR ODE INTEGRATION #####################################################

	// ********************************  8-th ORDER ADAPTIVE ADAMS-BASHFORTH-MOULTON STEPPER ********************************
	boost::numeric::odeint::adaptive_adams_bashforth_moulton<8UL, state_type, double, state_type, double,
															 boost::numeric::odeint::vector_space_algebra, boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer>
		abm8_stepper;

	// ********************************  4-th ORDER CLASSIC RUNGE-KUTTA STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta4_classic<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk4_stepper;

	// ********************************  5-th ORDER CASH-KARP STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rkk54_stepper;

	// ********************************  5-th ORDER DORMAND-PRINCE RUNGE-KUTTA ********************************
	// typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra> rkd5_stepper;

	// ********************************  BULIRSCH-STOER STEPPER ********************************
	// typedef boost::numeric::odeint::bulirsch_stoer<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> blstr_stepper;

	// ******************************** RUUNGE-KUTTA-FEHLBERG (RKF78) STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk78_stepper;

	// #################################################################################################################################################

	// start and end points, in terms of arc-length s, of each CTR segment and initial step-size for integration (ds)
	double s_start, s_end, ds;

	// instantiating the vector of initial conditions for solving the state equations (15 x 1)
	state_type y_0;

	/*
	 *****************************************************************************
	 * ========== initializing the initial conditions vector (15 x 1) ========== *
	 *****************************************************************************
	 */

	// 1st and 2nd element of y are the bending moments in the sheath along the x, y directions (in the local/body frame)
	y_0[0UL] = initGuess[0UL];
	y_0[1UL] = initGuess[1UL];
	// 3rd-5th elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
	y_0[2UL] = initGuess[2UL];
	y_0[3UL] = initGuess[3UL];
	y_0[4UL] = initGuess[4UL];
	// 6th-8th elements of y are twist angles, theta_i = [theta_1  theta_2  theta_3]
	y_0[5UL] = this->m_theta_0[0UL];
	y_0[6UL] = this->m_theta_0[1UL];
	y_0[7UL] = this->m_theta_0[2UL];
	// 9th-11th elements of y define the origin of the local frame propagates along the backbone at each arc-lengh 0 < s < L
	y_0[8UL] = 0.00;
	y_0[9UL] = 0.00;
	y_0[10UL] = 0.00;
	// 12th-15th elements of y are the (quaternion-orientations) of the local frame at each arc-length s
	y_0[11UL] = this->m_h_0[0UL];
	y_0[12UL] = this->m_h_0[1UL];
	y_0[13UL] = this->m_h_0[2UL];
	y_0[14UL] = this->m_h_0[3UL];

	// iterating through the tube segments comprising the CTR
	size_t len_seg = S.size() - 1UL;
	for (size_t seg = 0; seg < len_seg; ++seg)
	{
		// specifying the interval of integration (in terms of tube segment arc-lengths)
		s_start = S[seg];
		s_end = S[seg + 1UL];
		ds = (s_end - s_start) / 25.00; // 25 points per segment

		// passing the tube parameters in the segment to the state equation method
		this->m_stateEquations->setEquationParameters(blaze::column(U_x, seg), blaze::column(U_y, seg), blaze::column(EI, seg), blaze::column(GJ, seg), this->m_wf);

		// ##################################################### NUMERICAL INTEGRATION #####################################################
		// Employs the selected stepper (Numerical method) and integrates the system of ODEs along the segment considered
		boost::numeric::odeint::integrate_adaptive(abm8_stepper, *this->m_stateEquations, y_0, s_start, s_end, ds, *this->m_stateObserver);
	}

	//
	//	****************  #####  -------------- ___ DISTAL BOUNDARY CONDITIONS ___ --------------  #####  ****************
	//			1) internal moment at the tip of tube 1 must equal the external moment applied
	//			   Namely, mb_xy(L1) - (R1'L)_xy = 0
	//
	//			2) at the distal ends of the remaining tubes, the axial component of the internal moments must equal zero
	//			   Namely, ui_z(Li) - ui_z*(Li) = 0

	blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1;

	// grabbing the orientation at the distal end
	mathOp::getSO3(blaze::subvector<11UL, 4UL>(y_0), R1);

	blaze::StaticVector<double, 3UL> distalMoment = blaze::trans(R1) * m_wm;

	// Residue vector due to infringment of the distal boundary conditions || Residue = [ mb_x - Wm_x, mb_y - Wm_y, u1_z, u2_z, u3_z ]
	blaze::StaticVector<double, 5UL> Residue = {y_0[0UL] - distalMoment[0UL],
												y_0[1UL] - distalMoment[1UL],
												GJ(0UL, 0UL) * y_0[2UL] - blaze::trans(m_e3) * distalMoment,
												0.00,
												0.00};

	// lambda function that finds the u_z curvatures at the distal ends of tubes 2 and 3
	auto computeResidue = [&](double distalEnd, size_t index) -> void
	{
		// must use some tolerance when comparing floating points
		auto itt = std::lower_bound(this->m_s.begin(), this->m_s.end(), distalEnd - 1.00E-7); // finds where tube ends (with a 0.0001mm tolerance)

		auto id = std::distance(this->m_s.begin(), itt);
		// ui_z at the distal end of the i-th tube
		Residue[2UL + index] = this->m_y[id][2UL + index];
	};

	// Computing the Residues associated to the twist curvatures (rate of twist) at the distal ends of tubes 2 and 3
	blaze::StaticVector<double, 3UL> distEnd(this->m_segment->getDistalEnds()); // arc-lengths at which the distal ends of the tubes are currently

	computeResidue(distEnd[1UL], 1UL);
	computeResidue(distEnd[2UL], 2UL);

	return Residue;
}

// function that computes the finite-differences Jacobian for solving the BVP
blaze::StaticMatrix<double, 5UL, 5UL> CTR::jac_BVP(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 5UL> &residue)
{
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;

	blaze::StaticVector<double, 5UL> initGuessPerturbed(initGuess), residuePerturbed, scaled(initGuess);
	double incr_scale = 1.00E-7, incr_floor = 1.00E-9; // 1.00E-7, incr_floor = 1.00E-9 ==>> SEEM TO BE OPTIMAL;

	scaled *= incr_scale;
	scaled = blaze::generate(5UL, [&](size_t idx)
							 { return (std::fabs(scaled[idx]) > incr_floor) ? scaled[idx] : incr_floor; });

	for (size_t iter = 0UL; iter < 5UL; ++iter)
	{
		initGuessPerturbed[iter] += scaled[iter];
		// perturbed residue
		residuePerturbed = this->ODESolver(initGuessPerturbed);
		// building the finite-differences Residue jacobian
		blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / scaled[iter];
		// restoring the original value of the array
		initGuessPerturbed[iter] = initGuess[iter];
	}

	return jac_bvp;
}

// function that computes the finite-differences Jacobian wrt actuation inputs
blaze::StaticMatrix<double, 3UL, 6UL> CTR::jacobian(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 3UL> &tipPos)
{
	blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> jac;
	blaze::StaticVector<double, 6UL> q_Original(this->m_q), q_Perturbed(this->m_q), q_Scaled(this->m_q);
	double incr_scale = 1.00E-3, incr_floor = 5.00E-4;

	q_Scaled *= incr_scale;
	q_Scaled = blaze::generate(6UL, [&](size_t idx)
							   { return (std::fabs(q_Scaled[idx]) > incr_floor) ? q_Scaled[idx] : incr_floor; });

	for (size_t iter = 0UL; iter <= 5UL; ++iter)
	{
		q_Perturbed[iter] += q_Scaled[iter];
		this->setConfiguration(q_Perturbed);

		// recalculates the CTR transition points and segments for beta actuation
		if (iter <= 3UL)
			m_segment->recalculateSegments(this->m_Tubes, this->m_beta);

		this->ODESolver(initGuess);

		// computing the tip position of the perturbed CTR
		blaze::column(jac, iter) = (this->getTipPos() - tipPos) / q_Scaled[iter];
		// restoring the original CTR joint values
		q_Perturbed[iter] = q_Original[iter];
	}

	// sets the joint values to their original, unperturbed configuration values
	this->setConfiguration(q_Original);
	// this->ODESolver(initGuess);

	return jac;
}

// function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
bool CTR::PowellDogLeg(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	size_t k = 0UL;
	const size_t k_max = 300UL;
	double alpha, beta, delta, eps1, eps2, rho, c;
	blaze::StaticVector<double, 5UL> g, f, f_new, x_new, h_sd, h_gn, h_dl;
	blaze::StaticMatrix<double, 5UL, 5UL> J;

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	// initializing parameters
	delta = 1.00;
	eps1 = eps2 = 1.00e-22;

	f = this->ODESolver(initGuess);
	J = this->jac_BVP(initGuess, f);
	g = blaze::trans(J) * f;

	// checking if the initial guess satisfies the BVP without the need of any further refinement
	found = ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1)) ? true : false;

	while (!found && (k < k_max))
	{
		k++;

		alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
		h_sd = -alpha * g;			 // steepest descend (this is a direction, not a step!)
		h_gn = -mathOp::pInv(J) * f; // Gauss-Newton step (Least Square solution)

		// two candidates for the step to take from this point, a = alpha*h_sd & b = h_gn

		// computing the dog leg direction
		if (blaze::norm(h_gn) <= delta)
			h_dl = h_gn;
		else
		{
			if (blaze::norm(h_sd) >= delta)
				h_dl = delta * blaze::normalize(h_sd);
			else
			{
				c = blaze::trans(h_sd) * (h_gn - h_sd);

				if (c <= 0.00)
				{
					beta = (-c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) / blaze::sqrNorm(h_gn - h_sd);
				}
				else
				{
					beta = (delta * delta - blaze::sqrNorm(h_sd)) / (c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
				}

				h_dl = h_sd + beta * (h_gn - h_sd); // Dog Leg step
			}
		}

		if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
			found = true;
		else
		{
			x_new = initGuess + h_dl;

			f_new = this->ODESolver(x_new);
			rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.50 * blaze::trans(h_dl) * ((delta * h_dl) - g));

			if (rho > 0.00)
			{
				initGuess = std::move(x_new);
				f = std::move(f_new);
				J = this->jac_BVP(initGuess, f);
				g = blaze::trans(J) * f;

				if ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1))
					found = true;
			}

			if (rho > 0.75)
				delta = std::max(delta, 3.00 * blaze::norm(h_dl));
			else
			{
				if (rho < 0.25)
					delta *= 0.50;
			}

			if (delta < eps2 * (blaze::norm(initGuess) + eps2))
				found = true;
		}
	}

	return found;
}

// function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
bool CTR::Levenberg_Marquardt(blaze::StaticVector<double, 5UL> &initGuess)
{
	size_t k = 0UL;
	const size_t k_max = 300UL;
	blaze::StaticVector<double, 5UL> h, g, f, f_new;
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> J, A;
	blaze::IdentityMatrix<double> I(5UL);
	double rho, nu = 2.00, mu, tau = 1.00e-3, e1 = 1.00e-18, e2 = 1.00e-25;
	bool found;

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	// computing the residue and residue Jacobian associated to initGuess
	f = this->ODESolver(initGuess);
	J = this->jac_BVP(initGuess, f);
	A = blaze::trans(J) * J;
	g = blaze::trans(J) * f;
	found = (blaze::linfNorm(g) <= e1) ? true : false;
	mu = tau * blaze::max(blaze::diagonal(A));

	// starting the iterative minimization loop
	while ((!found) && (k < k_max))
	{
		k++;
		blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

		f_new = this->ODESolver(initGuess + h);
		rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.50 * blaze::trans(h) * ((mu * h) - g));

		if (rho > 0.00)
		{
			// accept the decrease in the function
			initGuess += h;
			// computing the residue Jacobian at the new initial guess
			J = this->jac_BVP(initGuess, f_new);
			A = blaze::trans(J) * J;
			f = std::move(f_new);
			g = blaze::trans(J) * f;
			found = (blaze::linfNorm(g) <= e1) ? true : false;
			mu = mu * std::max(0.33333333, 1.00 - blaze::pow(2.00 * rho - 1.00, 3.00));
			nu = 2.00;
		}
		else
		{
			mu = mu * nu;
			nu = 2.00 * nu;
		}

		// checking if the tolerance has been satisfied
		if (blaze::linfNorm(f) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Broyden (Nonlinear root-finding method for solving the BVP)
bool CTR::Broyden(blaze::StaticVector<double, 5UL> &initGuess)
{
	// found: returns true (false) when the root-finding method converges (does not converge) within k_max iterations
	bool found;

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = u3_z(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> JacInv, JacInvNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacInvNew = JacInv = mathOp::pInv(this->jac_BVP(X, F));

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		JacInv = std::move(JacInvNew);
		if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.00))
			JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) * blaze::trans(deltaX) * JacInv;
		else
			JacInvNew = JacInv;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - JacInv * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			X *= 0.75;
			readjustInitialGuesses(X);

			F = this->ODESolver(X);
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP(X, F));
			Xold = std::move(X);
			X = Xold - JacInv * F;
		}

		if (k % 10 == 0.00)
		{
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP(X, F));
			X = Xold - JacInv * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements Broyden's Nonlinear root-finding method for solving the BVP
bool CTR::Broyden_II(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> Jac, JacNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacNew = this->jac_BVP(initGuess, F);

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		Jac = std::move(JacNew);
		if (blaze::sqrNorm(deltaX) > 0.00)
			JacNew = Jac + blaze::sqrNorm(deltaX) * (deltaF - (Jac * deltaX)) * blaze::trans(deltaX);
		else
			JacNew = Jac;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - mathOp::pInv(Jac) * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			X *= 0.75;
			readjustInitialGuesses(X);

			F = this->ODESolver(X);
			JacNew = this->jac_BVP(X, F);
			Xold = std::move(X);
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (k % 10 == 0.00)
		{
			JacNew = this->jac_BVP(X, F);
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	// setting up and starting my handmade Newton-Raphson method
	blaze::StaticVector<double, 5UL> Residue, Residue_new, d_Residue, int_Residue, dGuess; // staticVectors are automatically initialized to 0

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = u3_z(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	// Residue of the unperturbed initial guess for the CTR
	Residue = this->ODESolver(initGuess);

	found = (blaze::linfNorm(Residue) <= this->m_accuracy) ? true : false;

	//  Weighing matrices for adjusting the initial guess iteratively (Implementing a PD regulator)
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 5UL, 5UL, blaze::rowMajor>> Kp, Ki, Kd;
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;
	blaze::diagonal(Kp) = 0.450; // 0.45 | 0.6  | 0.3
	blaze::diagonal(Ki) = 0.005;
	blaze::diagonal(Kd) = 0.002; // 3e-3 | 5e-3 | 2e-3

	size_t k = 0UL;
	const size_t k_max = 300UL;

	// starting iterations for adjusting the initial guess "u_guess ~ initGuess"
	while (!found && (k < k_max))
	{
		k++;
		jac_bvp = this->jac_BVP(initGuess, Residue);
		// error equation(globally asymptotically stable)
		dGuess = mathOp::pInv(jac_bvp) * (Kp * Residue + Ki * int_Residue + Kd * d_Residue);
		// updating the initial guess(weighted negative gradient of the cost function)
		initGuess -= dGuess;

		readjustInitialGuesses(initGuess);

		// computing the new cost associated to the newly readjusted initial guess
		Residue_new = this->ODESolver(initGuess);

		// cost variation due to initial guess refinement
		d_Residue = Residue_new - Residue;
		// integral of the residue
		int_Residue += Residue_new;
		// updating the cost
		Residue = std::move(Residue_new);

		if (blaze::linfNorm(Residue) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Modified_Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	/*
		Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition by Josef Stoer & Roland Bulirsch
	*/

	// zeroes |mb_x(0)|, |mb_y(0)|and limits the values of |u1_z(0)|, |u2_z(0)| and |u3_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses) -> void
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL] = 0.00;
	};

	readjustInitialGuesses(initGuess);

	bool found;
	// computes the residue associated to the initial guess
	blaze::StaticVector<double, 5UL> f(this->ODESolver(initGuess)), d;
	blaze::StaticVector<double, 5UL, blaze::rowVector> Dh;
	blaze::StaticMatrix<double, 5UL, 5UL> D, D_inv;
	double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
	size_t j = 0UL, k = 0UL;
	const size_t k_max = 300UL;
	std::vector<double> h_k; // vector to store all h_k's
	h_k.reserve(k_max);

	found = (blaze::linfNorm(f) <= this->m_accuracy) ? true : false;

	auto setupMethod = [&]() -> void
	{
		// then recomputes the residue
		f = this->ODESolver(initGuess);

		// computing the residue Jacobian
		D = this->jac_BVP(initGuess, f);

		// verifies NAN in the Jacobian and refines initial guess if necessary
		while (!blaze::isfinite(D))
		{
			initGuess *= 0.75;
			readjustInitialGuesses(initGuess);
			// then recomputes the residue
			f = this->ODESolver(initGuess);
			// computing the residue Jacobian
			D = this->jac_BVP(initGuess, f);
		}

		// pseudo-inverse of the residue Jacobian
		D_inv = mathOp::pInv(D);

		// search direction (directional derivative)
		d = D_inv * f;
		gamma = 1.00 / (blaze::norm(D_inv) * blaze::norm(D)); // gamma := 1/cond(Df)
		h_0 = blaze::sqrNorm(f);							  // h := f'f
		// Dh := D(f'f) = 2f'Df
		Dh = 2.00 * blaze::trans(f) * D;
		d_norm = blaze::norm(d);
		Dh_norm = blaze::norm(Dh);
	};

	while (!found && (k < k_max))
	{
		k++;
		setupMethod();

		while (true)
		{
			f = this->ODESolver(initGuess - blaze::pow(0.50, j) * d);
			// std::cout << "Modified_Newton_Raphson -- j = : " << j << " | residue = " << blaze::trans(f);
			while (!blaze::isfinite(f))
			{
				j++;
				f = this->ODESolver(initGuess - blaze::pow(0.50, j) * d);
				if (j > 10UL)
				{
					initGuess *= 0.75;
					readjustInitialGuesses(initGuess);
					setupMethod();
					j = 0UL;
				}
			}
			h = blaze::sqrNorm(f);
			improvementFactor = blaze::pow(0.50, j) * 0.25 * gamma * d_norm * Dh_norm;
			// storig the value of h_k to determine step size posteriorly
			h_k.push_back(h);

			if (h <= (h_0 - improvementFactor))
				break;
			else
				j++;
		}

		// retrieving the minimum h_k ==> h_k is monotonically decreasing (just grab its last element)
		lambda = blaze::pow(0.50, h_k.size() - 1UL);
		initGuess -= lambda * d;
		h_k.clear();

		// resets the exponent variable j
		j = 0UL;

		// checking the terminating condition
		if (blaze::linfNorm(f) <= this->m_accuracy)
		{
			return true;
		}
	}

	if (!found)
	{
		readjustInitialGuesses(initGuess);
		initGuess *= 0.75;
		found = this->PowellDogLeg(initGuess);

		if (!found)
		{
			readjustInitialGuesses(initGuess);
			initGuess *= 0.75;
			found = this->Levenberg_Marquardt(initGuess);
		}
	}

	return found;
}

// function that implements the CTR actuation for any inputs joint values q
bool CTR::actuate_CTR(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL> &q_input)
{
	// boolean flag for indicating convergence (1: zero found | 0: zero not found)
	bool found = false;

	// updating the CTR joints for desired input values
	this->setConfiguration(q_input);

	// recalculates the CTR transition points and segments
	m_segment->recalculateSegments(this->m_Tubes, this->m_beta);

	// initial guess for proximal boundary condition--[u_x(0) u_y(0) u1_z(0) u2_z(0) u3_z(0)]
	switch (this->m_method)
	{
	case mathOp::rootFindingMethod::NEWTON_RAPHSON:
		found = this->Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::LEVENBERG_MARQUARDT:
		found = this->Levenberg_Marquardt(initGuess);
		break;
	case mathOp::rootFindingMethod::POWELL_DOG_LEG:
		found = this->PowellDogLeg(initGuess);
		break;
	case mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON:
		found = this->Modified_Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN:
		found = this->Broyden(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN_II:
		found = this->Broyden_II(initGuess);
		break;
	}

	return found;
}

// function that implements the position control ==> returns timeout [bool]
bool CTR::posCTRL(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 3UL> &target, const double posTol)
{
	double minError = 1.00E3;										 // minimum distance to target
	bool status;													 // status = TRUE (FALSE) indicates convergence (lack thereof)
	blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> J;	 // Jacobian matrix
	blaze::StaticMatrix<double, 6UL, 3UL, blaze::columnMajor> J_inv; // Jacobian pseudoinverse
	blaze::IdentityMatrix<double, blaze::rowMajor> I(6UL);			 // 6 x 6 Identity matrix

	// zeroes |mb_x(0)|, |mb_y(0)| and |u3_z(0)| and limits the values of |u1_z(0)| and |u2_z(0)| to avoid numerical instability and lack of convergence
	auto readjustInitialGuesses = [](blaze::StaticVector<double, 5UL> &initial_guesses)
	{
		blaze::subvector<2UL, 3UL>(initial_guesses) = blaze::map(blaze::subvector<2UL, 3UL>(initial_guesses), [](double d)
																 { return (!blaze::isfinite(d)) ? 0.00 : blaze::sign(d) * std::min(blaze::abs(d), 50.00); });
		// mb_x(0) = mb_y(0) = 0.00;
		initial_guesses[0UL] = initial_guesses[1UL];
	};

	readjustInitialGuesses(initGuess);

	// proportional, derivative, and integral gains for position control
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>> Kp, Kd, Ki;
	blaze::diagonal(Kp) = 1.000; // 0.7500;
	blaze::diagonal(Ki) = 0.050; // 0.0005;
	blaze::diagonal(Kd) = 0.001; // 0.0001;	

	// Capturing the CTR's current joint configuration
	blaze::StaticVector<double, 6UL> dqdt, q_min(this->m_q), q(this->m_q);
	// Capturing the proximal BC for the minimum distance
	blaze::StaticVector<double, 5UL> initGuessMin(initGuess);
	// Calculate the CTR Jacobian in the present configuration and retrieves convergence status
	status = this->actuate_CTR(initGuess, q);

	// failure of convergence of the nonlinear root-finders
	if (!status)
		return status;

	blaze::StaticVector<double, 3UL> x_CTR = this->getTipPos();
    blaze::StaticVector<double, 3UL> tipError = target - x_CTR;
    blaze::StaticVector<double, 3UL> last_tipError = tipError;
    blaze::StaticVector<double, 3UL> d_tipError, int_tipError;

	// Euclidean distance to target
	double dist2Tgt = blaze::norm(tipError);

	if (dist2Tgt < minError)
	{
		minError = dist2Tgt;
		q_min = q;

		if (dist2Tgt <= posTol)
			return status;
	}

	// function to implement actuators sigularity avoidance
	blaze::StaticVector<double, 6UL> f;
	blaze::StaticVector<double, 3UL> f1 = blaze::subvector<0UL, 3UL>(f);

	// clearance between linear actuators
	double Clr = 5.00E-3, deltaBar = 0.00;

	// lengths of straight sections of the CTR tubes
	blaze::StaticVector<double, 3UL> ls, L;
	L[0UL] = this->m_Tubes[0UL]->getTubeLength();
	L[1UL] = this->m_Tubes[1UL]->getTubeLength();
	L[2UL] = this->m_Tubes[2UL]->getTubeLength();

	ls[0UL] = this->m_Tubes[0UL]->getStraightLen();
	ls[1UL] = this->m_Tubes[1UL]->getStraightLen();
	ls[2UL] = this->m_Tubes[2UL]->getStraightLen();
	// lower and upper bounds on prismatic joint limits
	blaze::StaticVector<double, 3UL> betaMax, betaMin;

	// iterations counter
	size_t N_itr = 0UL;
	// maximum admissible number of iterations in the position control loop
	const size_t maxIter = 500UL;
	// parameters for local optimization (joint limits avoidance)
	double ke = 2.00;

	// position control loop
	while ((dist2Tgt > posTol) && (N_itr < maxIter))
	{
		// incrementing the number of iterations
		N_itr++;

		// compute the Jacobian in the present configuration
		J = this->jacobian(initGuess, x_CTR);

		while (!blaze::isfinite(J))
		{
			initGuess *= 0.750;
			readjustInitialGuesses(initGuess);
			this->ODESolver(initGuess);
			x_CTR = this->getTipPos();
			J = this->jacobian(initGuess, x_CTR);
		}

		// Pseudo-inverse of Jacobian for resolving CTR joint motion rates
		J_inv = mathOp::pInv(J);

		// Nullspace control (collision and actuation limits)
		betaMin[0UL] = std::max({-ls[0UL] + deltaBar, L[1UL] + this->m_beta[1UL] - L[0UL], L[2UL] + this->m_beta[2UL] - L[0UL]});
		betaMin[1UL] = std::max({-ls[1UL] + deltaBar, this->m_beta[0UL] + Clr, L[2UL] + this->m_beta[2UL] - L[1UL]});
		betaMin[2UL] = std::max(-ls[2UL] + deltaBar, this->m_beta[1UL] + Clr);

		betaMax[0UL] = m_beta[1UL] - Clr;
		betaMax[1UL] = std::min(this->m_beta[2UL] - Clr, L[0UL] + this->m_beta[0UL] - L[1UL]);
		betaMax[2UL] = std::min({-deltaBar, L[1UL] + this->m_beta[1UL] - L[2UL], L[0UL] + this->m_beta[0UL] - L[2UL]});

		// penalty function for local optimization (actuator collision avoidance)
		f1 = blaze::pow(blaze::abs((betaMax + betaMin - 2.00 * this->m_beta) / (betaMax - betaMin + 1.00E-10)), ke) * blaze::sign(this->m_beta - (betaMax + betaMin) * 0.50);

		// resolved rates -- Nullspacec local optimization (joint limit avoidance)
		dqdt = J_inv * (Kp * tipError + Kd * d_tipError + Ki * int_tipError) + (I - blaze::trans(J_inv * J)) * (-f);

		auto rescale_dqdt = [&]() -> void 
		{ // rescaling linear joint variables for limit avoidance
			for (size_t i = 0UL; i < 3UL; ++i)
			{
				if (this->m_beta[i] + dqdt[i] > betaMax[i])
				{
					dqdt[i] = (betaMax[i] - this->m_beta[i]) * 0.50;
					// std::cout << "Entered hard limit constraint 1! |" << std::endl;
				}

				if (this->m_beta[i] + dqdt[i] < betaMin[i])
				{
					dqdt[i] = (betaMin[i] - this->m_beta[i]) * 0.50;
					// std::cout << "| Entered hard limit constraint 2!" << std::endl;
				}
			}
		};

		rescale_dqdt();

		// updating the CTR joints->q: [beta, theta]
		q += dqdt;

		// wrapping the actuation angles to the [0.00,2Pi) interval
		blaze::subvector<3UL, 3UL>(q) = blaze::map(blaze::subvector<3UL, 3UL>(q), [](double theta)
												   { return mathOp::congruentAngle(theta); });

		// actuate the CTR to new configuration and retrieve execution timeout status
		status = this->actuate_CTR(initGuess, q);

		// interrupts the loop execution if actuation fails
		if (!status)
		{
			initGuess *= 0.75;
			readjustInitialGuesses(initGuess);
			status = this->actuate_CTR(initGuess, q);

			if (!status)
			{
				// std::cout << "Nonlinear root-finders failed to converge: " << __PRETTY_FUNCTION__ << std::endl;
				initGuess = initGuessMin;
				this->actuate_CTR(initGuess, q_min);
				// return status;
			}
		}

		// tip position as predicted by the model
		x_CTR = this->getTipPos();

		// current position error
		tipError = target - x_CTR;
		// integrating the position error
		int_tipError += tipError;
		// derivative of the position error
		d_tipError = tipError - last_tipError;
		// updating the last tip error variable
		last_tipError = tipError;

		dist2Tgt = blaze::norm(tipError);

		if (dist2Tgt < minError)
		{
			minError = dist2Tgt;
			q_min = q;
			initGuessMin = initGuess;
		}

		// stops the control loop when the position update becomes significantly small
		if (blaze::linfNorm(dqdt) <= 1.00E-8)
		{
			initGuess = initGuessMin;
			status = this->actuate_CTR(initGuess, q_min);

			// std::cout << "Exited out of position control loop due small incremental threshold!" << std::endl;
			return status;
		}
	}

	// Actuating the CTR to the configuration which yields the minimum position error
	initGuess = std::move(initGuessMin);
	status = this->actuate_CTR(initGuess, q_min);

	return status;
}

// function that returns the Vector of tubes comprising the CTR
std::array<std::shared_ptr<Tube>, 3UL> CTR::getTubes()
{
	return this->m_Tubes;
}

// function that returns the current linear joint values of the CTR
blaze::StaticVector<double, 3UL> CTR::getBeta()
{
	return this->m_beta;
}

// function that returns the current joint values of the CTR
blaze::StaticVector<double, 6UL> CTR::getConfiguration()
{
	return this->m_q;
}

// function that returns the position of the CTR tip
blaze::StaticVector<double, 3UL> CTR::getTipPos()
{
	blaze::StaticVector<double, 3UL> pos;
	if (!this->m_y.empty())
		pos = blaze::subvector<8UL, 3UL>(this->m_y.back());

	return pos;
}

// function that returns the arc-lenghts at each tube's distal end
blaze::StaticVector<double, 3UL> CTR::getDistalEnds()
{
	return this->m_segment->getDistalEnds();
}

// function that returns the individual tube shapes
std::tuple<blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>> CTR::getTubeShapes()
{
	blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor> Tb_1(3UL, this->m_y.size());
	blaze::StaticVector<double, 3UL> v;

	// arc-lengths at the distal ends of each tube
	blaze::StaticVector<double, 3UL> distal_idx, distalEnds(this->m_segment->getDistalEnds());

	for (size_t col = 0UL; col < Tb_1.columns(); ++col)
	{
		v = blaze::subvector<8UL, 3UL>(this->m_y[col]);
		blaze::column(Tb_1, col) = v * 1.00E3;
	}

	// lambda returns the array index at which the ith tube ends
	auto tubeEndIndex = [&](size_t tube_index) -> size_t
	{
		// find the index in the arc-length vector at which each tube ends
		std::vector<double>::iterator it = std::lower_bound(this->m_s.begin(), this->m_s.end(), distalEnds[tube_index] - 1.00E-7); // finds where tube ends (0.0001mm tolerance)

		return std::distance(this->m_s.begin(), it);
	};

	size_t distalIndex_Tb2, distalIndex_Tb3;

	// index at which tube 2 ends
	distalIndex_Tb2 = tubeEndIndex(1UL);
	// index at which tube 3 ends
	distalIndex_Tb3 = tubeEndIndex(2UL);

	// number of columns in the shape matrices: index + 1
	auto Tb_2 = blaze::submatrix(Tb_1, 0UL, 0UL, 3UL, distalIndex_Tb2 + 1);
	auto Tb_3 = blaze::submatrix(Tb_1, 0UL, 0UL, 3UL, distalIndex_Tb3 + 1);

	// returns the tuple containing the shape of the tubes
	return std::make_tuple(Tb_1, Tb_2, Tb_3);
}

// function that returns a vector with the CTR shape
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> CTR::getShape()
{
	std::vector<double> r_x, r_y, r_z;
	r_x.reserve(this->m_y.size());
	r_y.reserve(this->m_y.size());
	r_z.reserve(this->m_y.size());

	if (this->m_y.size() > 0UL)
	{
		for (auto &el : this->m_y)
		{
			r_x.push_back(el[8UL]);
			r_y.push_back(el[9UL]);
			r_z.push_back(el[10UL]);
		}
	}

	return std::make_tuple(std::move(r_x), std::move(r_y), std::move(r_z));
}

// setter method for setting the actuation joint values (without actuating the CTR) <--> used for computing the Jacobian
void CTR::setConfiguration(const blaze::StaticVector<double, 6UL> &q)
{
	this->m_q = q;
	this->m_beta = blaze::subvector<0UL, 3UL>(this->m_q);
}

// function that sets which method to use for solving the BVP
void CTR::setBVPMethod(mathOp::rootFindingMethod mthd)
{
	this->m_method = mthd;
}

void CTR::setDistalMoment(const blaze::StaticVector<double, 3UL> &moment)
{
	this->m_wm = moment;
}

void CTR::setDistalForce(const blaze::StaticVector<double, 3UL> &force)
{
	this->m_wf = force;
}