#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include "Tube.hpp"

class Segment
{
private:
	std::vector<double> m_S;										  // arc-length of each tube transition point
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_EI;  // tubes bending stiffness in x,y directions
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_GJ;  // tubes torsional stiffness
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_U_x; // tubes' precurvature in the x  direction
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_U_y; // tubes' precurvature in the y direction
	std::array<double, 3UL> m_len_curv;								  // arc-length at which precurvature starts in the tubes
	std::array<double, 3UL> m_dist_end;								  // arc-length of tubes' distal ends

public:
	//default class constructor
	Segment();

	// overloaded class constructor
	Segment(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, const blaze::StaticVector<double, 3UL> &beta);

	// copy constructor
	Segment(const Segment &rhs);

	// move constructor
	Segment(Segment &&rhs) noexcept;

	// Segment desctructor
	~Segment() = default;

	// copy assignment operator
	Segment &operator=(const Segment &rhs);

	// move assignment operator
	Segment &operator=(Segment &&rhs) noexcept;

	// implements a functor to overload the constructors signature and allow recalculation of the CTR segmentation
	void recalculateSegments(const std::array<std::shared_ptr<Tube>, 3> &Tb, const blaze::StaticVector<double, 3UL> &beta);

	// getter method for retrieving the transition points defining the boundaries of all CTR segments
	std::vector<double> get_S();

	// getter method for returning the distal ends of all CTR tubes
	blaze::StaticVector<double, 3UL> getDistalEnds();

	// getter method for retrieving the vectors of tube bending stiffness in all CTR segments
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_EI();

	// getter method for retrieving the vectors of tube torional stiffness in all CTR segments
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_GJ();

	// getter method for retrieving the vectors of tube precurvatures along X in all CTR segments
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_U_x();

	// getter method for retrieving the vectors of tube precurvatures along Y in all CTR segments
	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_U_y();

	// method for returning all parameters along all CTR segments
	std::tuple<blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
			   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
			   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
			   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
			   std::vector<double>>
	returnParameters();
};