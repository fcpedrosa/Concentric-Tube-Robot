#pragma once

// #define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>
#include <blaze/Math.h>

/**
 * @brief Represents a tube in the concentric arrangement comprising the CTR robot
 */
class Tube
{
private:
	/**
	 * @brief Pi constant divided by 64.
	 */
	static const double pi_64;

	/**
	 * @brief Pi constant divided by 32.
	 */
	static const double pi_32;

	/**
	 * @brief Outer diameter of the tube.
	 */
	double m_OD;

	/**
	 * @brief Inner diameter of the tube.
	 */
	double m_ID;

	/**
	 * @brief Young's modulus of the tube material.
	 */
	double m_E;

	/**
	 * @brief Second moment of area of the tube cross-section.
	 */
	double m_I;

	/**
	 * @brief Shear modulus of the tube material.
	 */
	double m_G;

	/**
	 * @brief Polar moment of inertia of the tube cross-section.
	 */
	double m_J;

	/**
	 * @brief Stiffness matrix of the tube.
	 */
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> m_K;

	/**
	 * @brief Length of the straight segment of the tube.
	 */
	double m_ls;

	/**
	 * @brief Length of the curved segment of the tube.
	 */
	double m_lc;

	/**
	 * @brief Pre-curvature vector of the tube.
	 */
	blaze::StaticVector<double, 3UL> m_u_ast;

public:
	/**
	 * @brief Implements the default constructor for the Tube class
	 */
	Tube();

	/**
	 * @brief Implements the overloaded constructor for the Tube class.
	 *
	 * @param OD The outer diameter of the tube.
	 * @param ID The inner diameter of the tube.
	 * @param E The Young's modulus of the tube material.
	 * @param G The Shear modulus of the tube material
	 * @param ls The length of the straight transmission of the tube
	 * @param lc The length of the curved segment of the tube
	 * @param u_ast A 3-dimensional static Blaze vector of the pre-curvature of the tube
	 */
	Tube(double OD, double ID, double E, double G, double ls, double lc, const blaze::StaticVector<double, 3UL> &u_ast);

	/**
	 * @brief Destroys the Tube object.
	 */
	~Tube() = default;

	/**
	 * @brief Implements the copy constructor for the Tube class.
	 *
	 * @param rhs The source Tube object to copy from.
	 */
	Tube(const Tube &rhs);

	/**
	 * @brief Implements the move constructor for the Tube class.
	 *
	 * @param rhs The source Tube object to move from.
	 */
	Tube(Tube &&rhs) noexcept;

	/**
	 * @brief Implements the copy assignment operator for the Tube class.
	 *
	 * @param rhs The source Tube object to copy from.
	 * @return A reference to the assigned Tube object.
	 */
	Tube &operator=(const Tube &rhs);

	/**
	 * @brief Implements the move assignment operator for the Tube class.
	 *
	 * @param rhs The source Tube object to move from.
	 * @return A reference to the assigned Tube object.
	 */
	Tube &operator=(Tube &&rhs) noexcept;

	/**
	 * @brief Implements getter method for retrieving the tube's kinematic parameters
	 *
	 * @return A tuple containing the parameters [OD, ID, E, G, ls, lc, ||u_ast||] .
	 */
	std::tuple<double, double, double, double, double, blaze::StaticVector<double, 3UL>> getTubeParameters();

	/**
	 * @brief Implements a setter method for setting the Young's modulus of the Tube object
	 *
	 * @param E The numerical value of the new Young's modulus of the Tube's material
	 */
	void setYoungModulus(double E);

	/**
	 * @brief Implements a setter method for setting the Shear modulus of the Tube object
	 *
	 * @param G The numerical value of the new Shear modulus of the Tube's material
	 */
	void setShearModulus(double G);

	/**
	 * @brief Implements a getter method for retrieving the tube length
	 *
	 * @return The overall tube length in meters (straight + curved sections)
	 */
	double getTubeLength();

	/**
	 * @brief Implements a getter method for retrieving the pre-curvature of the Tube object
	 *
	 * @return A 3-dimensional static Blaze vector with the pre-curvature of the tube
	 */
	blaze::StaticVector<double, 3UL> get_u_ast();

	/**
	 * @brief Implements a getter method for retrieving the pre-curvature of the Tube object along a specific direction
	 *
	 * @return A scalar with the pre-curvature of the tube along the 'x' or 'y' directions
	 */
	double get_u_ast(const size_t id);

	/**
	 * @brief Implements a setter method for updating the pre-curvature of the Tube object
	 *
	 * @param u_ast 3-dimensional static Blaze vector with the new pre-curvature of the tube
	 */
	void set_u_ast(const blaze::StaticVector<double, 3UL> &u_ast);

	/**
	 * @brief Implements a setter method for updating the pre-curvature of the Tube object along a specific direction
	 *
	 * @param id A size_t index of the direction 0: 'x', 1: 'y'
	 * @param u scalar with the pre-curvature of the tube along the 'x' or 'y' directions, as determined by id
	 */
	void set_u_ast(const size_t id, const double u);

	/**
	 * @brief Implements a getter method for retrieving the length of the straight tranmission of the Tube object
	 *
	 * @return A scalar with the length of the straight segment of the tube in meters
	 */
	double getStraightLen();

	/**
	 * @brief Implements a getter method for retrieving the length of the curved segment of the Tube object
	 *
	 * @return A scalar with the length of the curved segment of the tube in meters
	 */
	double getCurvLen();

	/**
	 * @brief Implements a setter method for updating the length of the straight tranmission of the Tube object
	 *
	 * @param ls scalar with the new length of the straight segment of the tube in meters
	 */
	void setStraightLen(double ls);

	/**
	 * @brief Implements a setter method for updating the length of the curved segment of the Tube object
	 *
	 * @param lc scalar with the length of the curved segment of the tube in meters
	 */
	void setCurvLen(double lc);

	/**
	 * @brief Implements a getter method for retrieving the bending stiffness matrix for the Tube object
	 *
	 * @return A 3x3 static diagonal Blaze matrix witht the bending stiffness of the Tube object
	 */
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> getK_Matrix();

	/**
	 * @brief Implements a getter method for retrieving the bending/torsional stiffness of the Tube along a specific direction
	 *
	 * @param i size_t index of the direction: 0: 'x', 1: 'y', 2: 'z'
	 */
	double getK(int i);

	/**
	 * @brief Implements a setter method for updating the bending & torional stiffness of the Tube
	 *
	 * @param EI The new bending stiffness along the 'x' and 'y' directions
	 * @param GJ The new torsional stiffness along the 'z' direction
	 */
	void setK(const double EI, const double GJ);

	/**
	 * @brief Implements a setter method for updating the bending stiffness of the Tube
	 *
	 * @param EI The new bending stiffness along the 'x' and 'y' directions
	 */
	void setBendingK(const double EI);
};