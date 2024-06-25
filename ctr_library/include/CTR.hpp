// ***************************************************************************************** //
// *  This file is part of NIH_CSTAR, a CTR kinematics library for Concentric Tube Robots  * //
// *																					   * //
// *  ----------- # Copyright (C) 2021 Filipe C. Pedrosa <fpedrosa@uwo.ca> # -----------   * //
// *																					   * //
// *  Project developed under the supervision of Prof Dr Rajni Patel <rvpatel@uwo.ca>	   * //
// *			  CSTAR (Canadian Surgical Technologies & Advanced Robotics)			   * //
// *				   Western University, London, ON, Canada							   * //
// ***************************************************************************************** //

#pragma once

#include "Tube.hpp"
#include "Segment.hpp"
#include "ODESystem.hpp"
#include "Observer.hpp"
#include "mathOperations.hpp"
#include <boost/numeric/odeint.hpp>
#include <memory>
#include <tuple>

/**
 * @brief Class implementing a three-tube Concentric Tube Robot (CTR).
 */
class CTR
{
public:
	/**
	 * @brief Removes the default CTR Class constructor
	 */
	CTR() = delete;

	/**
	 * @brief Implements the overloaded constructor for the CTR class.
	 *
	 * @param Tb A 3-dimensional std::array containing smart pointers to the three tube objects comprising the CTR assembly.
	 * @param q A 6-dimensional static Blaze vector containing the actuation input values for the CTR.
	 * @param Tol A scalar that dictates the prescribed accuracy (tolerance) for solving the associated Boundary Value Problem (BVP).
	 * @param method An enum class specifying the root-finding method. The possible choices are: NEWTON_RAPHSON, LEVENBERG_MARQUARDT, POWELL_DOG_LEG, MODIFIED_NEWTON_RAPHSON, BROYDEN, and BROYDEN_II.
	 */
	CTR(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, blaze::StaticVector<double, 6UL> &q, const double Tol, const mathOp::rootFindingMethod method);

	/**
	 * @brief Implements the copy constructor for the CTR class.
	 *
	 * @param rhs The source CTR object to copy from.
	 */
	CTR(const CTR &rhs);

	/**
	 * @brief Implements the move constructor for the CTR class.
	 *
	 * @param rhs The source CTR object to move from.
	 */
	CTR(CTR &&rhs) noexcept;

	/**
	 * @brief Destroys the CTR object.
	 */
	~CTR() = default;

	/**
	 * @brief Implements the copy assignment operator for the CTR class.
	 *
	 * @param rhs The source CTR object to copy from.
	 * @return A reference to the assigned CTR object.
	 */
	CTR &operator=(const CTR &rhs);

	/**
	 * @brief Implements the move assignment operator for the CTR class.
	 *
	 * @param rhs The source CTR object to move from.
	 * @return A reference to the assigned CTR object.
	 */
	CTR &operator=(CTR &&rhs) noexcept;

	/**
	 * @brief Resets the std::vector containing all state variables for the CTR at each arc-length
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 */
	void reset(const blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Solves the CTR model equations (system of ODEs) by numerical integration
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a 5-dimensional residue vector. The residue is the violation of the distal boundary condition
	 */
	blaze::StaticVector<double, 5UL> ODESolver(const blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Computes the finite-differences Jacobian for solving the associated Boundary Value Problem (BVP)
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @param residue is the residue (violation of the distal boundary condition) at the current iteration of the shooting problem.
	 * @return a 5x5 BVP Jacobian matrix that relates how the residue changes as a function of the initial guesses
	 */
	blaze::StaticMatrix<double, 5UL, 5UL> jac_BVP(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 5UL> &residue);

	/**
	 * @brief Computes the finite-differences Jacobian for the spatial velocities at the distal-end of the CTR.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @param tipPos is the current distal-end position at the current iteration of the inverse kinematics (position control problem).
	 * @return a 3x6 BVP Jacobian matrix that relates how the distal spatial velocities residue change as a function of the actuation inputs
	 */
	blaze::StaticMatrix<double, 3UL, 6UL> jacobian(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 3UL> &tipPos);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the Powell's Dog Leg Method to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool PowellDogLeg(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the Levenberg-Marquardt Method to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool Levenberg_Marquardt(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the Broyden's Method (Jacobian inverse is estimated) to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool Broyden(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the Broyden's Method (Jacobian inverse is estimated) to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool Broyden_II(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the Newton-Raphson's Method to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Nonlinear root-finding algorithm for solving the shooting method. It implements the modified Newton-Raphson's Method to find the zero of the residue function.
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @return a boolean flag. True: the zero (root) of the residue function has been found | False: the zero (root) of the residue function has not been found
	 */
	bool Modified_Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess);

	/**
	 * @brief Actuates the CTR robot to a configuration determined by the actuation inputs and the associated boundary conditions dictated by the initial guess for the BVP
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @param q_input is a 6-dimensional vector of the actuation units The first three entries are the linear joints, the last three are the revolute ones.
	 * @return a boolean flag. True: the associated BVP has been solved successfully | False: otherwise
	 */
	bool actuate_CTR(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL> &q_input);

	/**
	 * @brief Implements the differential inverse kinematics based on resolved rates motion of the actuators
	 *
	 * @param initGuess is the initial guess 5-dimensional vector for the boundary value problem (BVP).
	 * @param target is a 3-dimensional vector of the target position to where the end-effector (distal-end of the CTR) should be steered to.
	 * @param Tol is a scalar indicating the tolerance with which the position control problem should be solved.
	 * @return a boolean flag. True: the associated IK has been solved successfully within the prescribed tolerance | False: otherwise
	 */
	bool posCTRL(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 3UL> &target, const double Tol);

	/**
	 * @brief Returns an std container with smart pointers to the Tube objects comprising the CTR
	 *
	 * @return smart pointers to the tubes comprising the CTR assembly
	 */
	std::array<std::shared_ptr<Tube>, 3UL> getTubes();

	/**
	 * @brief Returns the current actuation values of the linear joints of the CTR
	 *
	 * @return a static 3-dimensional Blaze vector with the actuation values for the linear joints [beta_1, beta_2, beta_3]
	 */
	blaze::StaticVector<double, 3UL> getBeta();

	/**
	 * @brief Returns the current actuation values of all the linear and revolute joints of the CTR
	 *
	 * @return a static Blaze vector with the actuation values for all the joints [beta_1, beta_2, beta_3, alpha_1, alpha_2, alpha_3]
	 */
	blaze::StaticVector<double, 6UL> getConfiguration();

	/**
	 * @brief Returns the current end-effector (distal-end) position of the CTR in meters
	 *
	 * @return a static 3-dimensional Blaze vector with the position of the distal-end of the CTR
	 */
	blaze::StaticVector<double, 3UL> getTipPos();

	/**
	 * @brief Returns the arc-lenghts of each tube's distal end in meters
	 *
	 * @return a static 3-dimensional Blaze vector with the arc-length distal-end positions of each tube in the CTR assembly
	 */
	blaze::StaticVector<double, 3UL> getDistalEnds();

	/**
	 * @brief Implements a getter moethod for acquiring the 3-dimensional shape of each tube in the CTR assembly
	 *
	 * @return a set of three 3xN Blaze matrices with the shape of each tube in the CTR assembly. The first, second, and third rows of each matrix correspond to the x,y,z coordinates
	 */
	std::tuple<blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>> getTubeShapes();

	/**
	 * @brief Returns a tuple with three std vectors containing the shape of the innermost tube in the CTR assembly
	 *
	 * @return a tuple with the set of three std::vectors with the shape of the innermost tube in the CTR assembly. The first, second, and third entries in the tuple correspond to the x,y,z coordinates of the tube centerline
	 */
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> getShape();

	/**
	 * @brief Implements a setter method for setting the actuation joint values. It does not actuate the CTR (does not solve the FK of the CTR)
	 *
	 * @param q is a static 6-dimensional Blaze vector containing the actuation values for the CTR joints [beta_1, beta_2, beta_3, alpha_1, alpha_2, alpha_3]
	 */
	void setConfiguration(const blaze::StaticVector<double, 6UL> &q);

	/**
	 * @brief Implements a setter method for chosing wich nonlinear root finding algorithm to use when solving the Boundary Value Problem (BVP)
	 *
	 * @param An enum class specifying the root-finding method. The possible choices are: NEWTON_RAPHSON, LEVENBERG_MARQUARDT, POWELL_DOG_LEG, MODIFIED_NEWTON_RAPHSON, BROYDEN, and BROYDEN_II.
	 */
	void setBVPMethod(mathOp::rootFindingMethod mthd);

	/**
	 * @brief Implements a setter method for specifying the point moment acting at the distal-end of the CTR
	 *
	 * @param moment is a static 3-dimensional Blaze vector containing the moment acting at the distal-end of the CTR.
	 */
	void setDistalMoment(const blaze::StaticVector<double, 3UL> &moment);

	/**
	 * @brief Implements a setter method for specifying the point force acting at the distal-end of the CTR
	 *
	 * @param force is a static 3-dimensional Blaze vector containing the point force acting at the distal-end of the CTR.
	 */
	void setDistalForce(const blaze::StaticVector<double, 3UL> &force);

private:
	/**
	 * @brief Defines the accuracy of the numerical solution for the CTR Boundary Value Problem (BVP).
	 */
	double m_accuracy;

	/**
	 * @brief Available root-finding methods: Newton-Raphson, Levenberg-Marquardt, Powell (dog-leg), Broyden.
	 */
	mathOp::rootFindingMethod m_method;

	/**
	 * @brief Array of smart pointers to the Tube objects comprising the CTR.
	 */
	std::array<std::shared_ptr<Tube>, 3UL> m_Tubes;

	/**
	 * @brief Linear actuation values [beta_1, beta_2, beta_3].
	 */
	blaze::StaticVector<double, 3UL> m_beta;

	/**
	 * @brief Joint actuation values [alpha_1, alpha_2, alpha_3].
	 */
	blaze::StaticVector<double, 6UL> m_q;

	/**
	 * @brief Initial twist angle for all tubes at s = 0.
	 */
	blaze::StaticVector<double, 3UL> m_theta_0;

	/**
	 * @brief Initial orientation of the local frame at s = 0, represented as a quaternion.
	 */
	blaze::StaticVector<double, 4UL> m_h_0;

	/**
	 * @brief External point force at the distal-end of the CTR (force component of the external wrench).
	 */
	blaze::StaticVector<double, 3UL> m_wf;

	/**
	 * @brief External point moment at the distal-end of the CTR (moment component of the external wrench).
	 */
	blaze::StaticVector<double, 3UL> m_wm;

	/**
	 * @brief Third canonical basis of R^3.
	 */
	blaze::StaticVector<double, 3UL> m_e3;

	/**
	 * @brief Unique pointer to the segments between transition points along the CTR backbone.
	 */
	std::unique_ptr<Segment> m_segment;

	/**
	 * @brief Stores the CTR state vector at each arc-length (spatial integration step).
	 */
	std::vector<state_type> m_y;

	/**
	 * @brief Stores the discrete arc-length points (in meters) along the backbone.
	 */
	std::vector<double> m_s;

	/**
	 * @brief Unique pointer to the ODE system implementing the state differential equations for a three-tube CTR.
	 */
	std::unique_ptr<ODESystem> m_stateEquations;

	/**
	 * @brief Unique pointer to the state observer for Boost::odeInt.
	 */
	std::unique_ptr<Observer> m_stateObserver;
};