// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Segment.hpp"

//default class constructor
Segment::Segment()
{
	m_S.clear();
	m_EI.reset();
	m_GJ.reset();
	m_U_x.reset();
	m_U_y.reset();
	m_len_curv = {0.00, 0.00, 0.00};
	m_dist_end = {0.00, 0.00, 0.00};
}

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

Segment::~Segment()
{
	// nothing to be done
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
	// vector of overall length of each tube
	std::array<double, 3UL> tb_len;
	// clearing tube transition points
	this->m_S.clear();

	for (size_t i = 0UL; i < 3UL; ++i)
	{
		// ith tube overall length
		tb_len[i] = Tb[i]->getTubeLength();
		// arc-length at which the ith tube ends
		this->m_dist_end[i] = tb_len[i] + beta[i];
		// arc-length at which the curved segment of each tube starts
		this->m_len_curv[i] = this->m_dist_end[i] - Tb[i]->getCurvLen();
	}

	const double TOLERANCE = 1.00E-7;
	auto compare = [&](double a, double b) -> bool
	{
		return std::fabs(a - b) < TOLERANCE;
	};

	// listing all landmark points at tube segment transitions (s >= 0)
	this->m_S.push_back(0.00);
	// inserting m_len_curv array into landmarks
	this->m_S.insert(this->m_S.end(), this->m_len_curv.begin(), this->m_len_curv.end());
	// inserting m_dist_end array into landmarks
	this->m_S.insert(this->m_S.end(), this->m_dist_end.begin(), this->m_dist_end.end());
	// sorting the landmark points
	std::sort(this->m_S.begin(), this->m_S.end());
	// acquiring the unique landmark points only
	auto it = std::unique(this->m_S.begin(), this->m_S.end(), compare); // tolerance used to allow for small differences due to precision limitations
	
	// deleting any repeated elements
	this->m_S.erase(it, this->m_S.end());

	// total number of segments in the current CTR configuration (transition points - 1)
	size_t len = this->m_S.size() - 1UL;
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

	// determining the indexes correponding to tube transitions
	for (size_t i = 0UL; i < 3UL; ++i)
	{
		// transition points (straight -> curved sections)
		element = this->m_len_curv[i];
		// finds where curved section of the i-th starts
		it_b = std::lower_bound(this->m_S.begin(), this->m_S.end(), element - TOLERANCE); // tolerance used to allow for small differences due to precision limitations

		// transition points -> distal ends
		element = m_dist_end[i];
		// finds where i-th tube ends
		it_c = std::lower_bound(this->m_S.begin(), this->m_S.end(), element - TOLERANCE); // tolerance used to allow for small differences due to precision limitations

		b = std::distance(this->m_S.begin(), it_b);
		c = std::distance(this->m_S.begin(), it_c);

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
std::vector<double> Segment::get_S()
{
	return this->m_S;
}

// getter method for returning the distal ends of all CTR tubes
blaze::StaticVector<double, 3UL> Segment::getDistalEnds()
{
	blaze::StaticVector<double, 3UL> distalEnds = { this->m_dist_end[0UL], this->m_dist_end[1UL], this->m_dist_end[2UL] };
	return distalEnds;
}

// getter method for retrieving the vectors of tube bending stiffness in all CTR segments
blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> Segment::get_EI()
{
	return this->m_EI;
}

// getter method for retrieving the vectors of tube torional stiffness in all CTR segments
blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> Segment::get_GJ()
{
	return this->m_GJ;
}

// getter method for retrieving the vectors of tube precurvatures along X in all CTR segments
blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> Segment::get_U_x()
{
	return this->m_U_x;
}

// getter method for retrieving the vectors of tube precurvatures along Y in all CTR segments
blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> Segment::get_U_y()
{
	return this->m_U_y;
}

// method for returning all parameters along all CTR segments
std::tuple<blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
		   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
		   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
		   blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
		   std::vector<double>>
Segment::returnParameters()
{

	return std::make_tuple(this->m_EI, this->m_GJ, this->m_U_x, this->m_U_y, this->m_S);
}
