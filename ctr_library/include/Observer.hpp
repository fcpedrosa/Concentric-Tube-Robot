#pragma once

#include <blaze/Math.h>
#include <vector>

typedef blaze::StaticVector<double, 15UL> state_type;

class Observer
{
private:
	std::vector<state_type> &m_states; // vector with all CTR states along entire backbone
	std::vector<double> &m_arcLength;  // vector with all discretized arc-length points along the backbone (CTR shape)

public:
	// defult constructor
	Observer() = delete;

	// overloaded class constructor
	Observer(std::vector<state_type> &states, std::vector<double> &s) : m_states(states), m_arcLength(s) {}

	// observer functor for capturing data from Boost::odeInt integrator
	void operator()(const state_type &states, double s);
};