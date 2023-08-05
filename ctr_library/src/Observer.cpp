// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Observer.hpp"
#include <iostream>

// observer functor for capturing data from Boost::odeInt integrator
void Observer::operator()(const state_type &states, double s)
{
	m_states.push_back(states);
	m_arcLength.push_back(s);
}