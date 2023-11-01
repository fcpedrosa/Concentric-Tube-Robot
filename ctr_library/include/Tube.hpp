#pragma once

// #define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>
#include <blaze/Math.h>

class Tube
{
private:
	static const double pi_64;
	static const double pi_32;
	double m_OD;
	double m_ID;
	double m_E;
	double m_I;
	double m_G;
	double m_J;
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> m_K;
	double m_ls;
	double m_lc;
	blaze::StaticVector<double, 3UL> m_u_ast;

public:
	// default constructor
	Tube();

	// Tube constructor
	Tube(double OD, double ID, double E, double G, double ls, double lc, const blaze::StaticVector<double, 3UL> &u_ast);

	// Tube desctructor
	~Tube() = default;

	// copy constructor
	Tube(const Tube &rhs);

	// move constructor
	Tube(Tube &&rhs) noexcept;

	// Copy assignment operator
	Tube &operator=(const Tube &rhs);

	// move assignment operator
	Tube &operator=(Tube &&rhs) noexcept;

	// get method for retrieving the tube parameters
	std::tuple<double, double, double, double, double, blaze::StaticVector<double, 3UL>> getTubeParameters();

	// set method for updating the Young's modulus
	void setYoungModulus(double E);

	// set method for updating the Shear modulus
	void setShearModulus(double G);

	// get method for retrieving the tube overall length
	double getTubeLength();

	// get method for retrieving the tube precurvature vector
	blaze::StaticVector<double, 3UL> get_u_ast();

	// get method for retrieving the "scalar" curvature along x or y directions
	double get_u_ast(const size_t id);

	// set method for updating the tube precurvature vector
	void set_u_ast(const blaze::StaticVector<double, 3UL> &u_ast);

	// set method for updadting the tube's precurvature along x or y directions
	void set_u_ast(const size_t id, const double u);

	// get method for retrieving the length of the straight section
	double getStraightLen();

	// get method for retrieving the length of the tube curved section
	double getCurvLen();

	// set method for updating the length of the straight section
	void setStraightLen(double ls);

	// set method for updating the legnth of the curved section
	void setCurvLen(double lc);

	// get method for retrieving the stiffness matrix
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> getK_Matrix();

	// get method for retrieving the ith entry of the main diagonal
	double getK(int i);

	// method for setting the bending & torsional stiffness
	void setK(const double EI, const double GJ);

	// method for setting the bending stiffness in the stiffness matrix
	void setBendingK(const double EI);
};