// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Segment.hpp"
#include <algorithm>

// overloaded constructor
Segment::Segment(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, const blaze::StaticVector<double, 3UL> &beta)
{
	this->recalculateSegments(Tb, beta);
}

// copy constructor
Segment::Segment(const Segment &rhs) : m_S(rhs.m_S), m_EI(rhs.m_EI), m_GJ(rhs.m_GJ), m_U_x(rhs.m_U_x),
									   m_U_y(rhs.m_U_y), m_len_curv(rhs.m_len_curv), m_dist_end(rhs.m_dist_end) {}

// move constructor
Segment::Segment(Segment &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_S = std::move(rhs.m_S);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_U_x = std::move(rhs.m_U_x);
		this->m_U_y = std::move(rhs.m_U_y);
		this->m_len_curv = std::move(rhs.m_len_curv);
		this->m_dist_end = std::move(rhs.m_dist_end);
	}
}

// copy assignment operator
Segment &Segment::operator=(const Segment &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_S = rhs.m_S;
		this->m_EI = rhs.m_EI;
		this->m_GJ = rhs.m_GJ;
		this->m_U_x = rhs.m_U_x;
		this->m_U_y = rhs.m_U_y;
		this->m_len_curv = rhs.m_len_curv;
		this->m_dist_end = rhs.m_dist_end;
	}
	return *this;
}

// move assignment operator
Segment &Segment::operator=(Segment &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_S = std::move(rhs.m_S);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_U_x = std::move(rhs.m_U_x);
		this->m_U_y = std::move(rhs.m_U_y);
		this->m_len_curv = std::move(rhs.m_len_curv);
		this->m_dist_end = std::move(rhs.m_dist_end);
	}
	return *this;
}

// implements a functor to overload the constructors signature and allow recalculation of the CTR segmentation
void Segment::recalculateSegments(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, const blaze::StaticVector<double, 3UL> &beta)
{
	// Vector of overall length of each tube
    const blaze::StaticVector<double, 3UL> tb_len = {
        Tb[0UL]->getTubeLength(),
        Tb[1UL]->getTubeLength(),
        Tb[2UL]->getTubeLength()
    };

	// clearing tube transition points
	this->m_S.clear();
	// Reserving a fixed capacity to avoid reallocations
	m_S.reserve(10UL);
	// listing all landmark points at tube segment transitions (s >= 0)
	this->m_S.emplace_back(0.00);

	// Arc-length at which each tube ends and the curved segment of each tube starts
    for (size_t i = 0; i < 3UL; ++i) 
	{
        this->m_dist_end[i] = tb_len[i] + beta[i];
        this->m_len_curv[i] = this->m_dist_end[i] - Tb[i]->getCurvLen();

		// Inserting segment transition points and tube end points
		this->m_S.emplace_back(this->m_len_curv[i]);
        this->m_S.emplace_back(this->m_dist_end[i]);
    }

	static constexpr double TOLERANCE = 1.00E-7;
	auto compare = [&](double a, double b) -> bool
	{
		return std::fabs(a - b) < TOLERANCE;
	};

	// Sorting and removing duplicates
	std::sort(this->m_S.begin(), this->m_S.end());
	auto it = std::unique(this->m_S.begin(), this->m_S.end(), compare);
	this->m_S.erase(it, this->m_S.end());

	// total number of segments in the current CTR configuration (transition points - 1)
	const size_t len = this->m_S.size() - 1UL;

	// Alocatting memory space for the output matrices
	this->m_EI.resize(3UL, len, false);
	this->m_GJ.resize(3UL, len, false);
	this->m_U_x.resize(3UL, len, false);
	this->m_U_y.resize(3UL, len, false);

	// filling matrices with defaul zero values
	this->m_EI = this->m_GJ = this->m_U_x = this->m_U_y = 0.00;

	size_t b, c, span;
	double element;
	std::vector<double>::iterator it_b, it_c;
	auto sBegin = this->m_S.begin();

	// determining the indexes correponding to tube transitions
	for (size_t i = 0UL; i < 3UL; ++i)
	{
		// transition points (straight -> curved sections)
		element = this->m_len_curv[i];
		// finds where curved section of the i-th starts
		it_b = std::lower_bound(sBegin, this->m_S.end(), element - TOLERANCE); // tolerance used to allow for small differences due to precision limitations

		// transition points -> distal ends
		element = m_dist_end[i];
		// finds where i-th tube ends
		it_c = std::lower_bound(sBegin, this->m_S.end(), element - TOLERANCE); // tolerance used to allow for small differences due to precision limitations

		b = std::distance(sBegin, it_b);
		c = std::distance(sBegin, it_c);

		// populating the output matrices accordingly
		blaze::submatrix(this->m_EI, i, 0UL, 1UL, c) = Tb[i]->getK(1); // bending stiffness
		blaze::submatrix(this->m_GJ, i, 0UL, 1UL, c) = Tb[i]->getK(3); // torsional stiffness

		// only load the curvatures when the tubes have a curved section
		if (b != c)
		{
			span = c - b;
			// loading the precurvatures along the x and y directionss
			blaze::submatrix(this->m_U_x, i, b, 1UL, span) = Tb[i]->get_u_ast(1); // precurvature along the x direction
			blaze::submatrix(this->m_U_y, i, b, 1UL, span) = Tb[i]->get_u_ast(2); // precurvature along the y direction
		}
	}
}

// getter method for retrieving the transition points defining the boundaries of all CTR segments
const std::vector<double> & Segment::get_S() const
{
	return this->m_S;
}

// getter method for returning the distal ends of all CTR tubes
const blaze::StaticVector<double, 3UL>& Segment::getDistalEnds() const
{
	return this->m_dist_end;
}

// getter method for retrieving the vectors of tube bending stiffness in all CTR segments
const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &Segment::get_EI() const
{
	return this->m_EI;
}

// getter method for retrieving the vectors of tube torional stiffness in all CTR segments
const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &Segment::get_GJ() const
{
	return this->m_GJ;
}

// getter method for retrieving the vectors of tube precurvatures along X in all CTR segments
const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &Segment::get_U_x() const
{
	return this->m_U_x;
}

// getter method for retrieving the vectors of tube precurvatures along Y in all CTR segments
const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &Segment::get_U_y() const
{
	return this->m_U_y;
}

// method for returning all parameters along all CTR segments
std::tuple<const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &,
		   const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &,
		   const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &,
		   const blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> &,
		   const std::vector<double> &>
Segment::returnParameters() const
{

	return std::tie(this->m_EI, this->m_GJ, this->m_U_x, this->m_U_y, this->m_S);
}