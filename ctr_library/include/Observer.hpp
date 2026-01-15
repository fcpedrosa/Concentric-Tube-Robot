#pragma once

#include <blaze/Math.h>
#include <vector>

typedef blaze::StaticVector<double, 15UL> state_type;

/**
 * @file Observer.hpp
 * @brief ODE integration observer utilities.
 * @ingroup ctr_utils
 * @details Captures states and arc-length samples during integration.
 */

/**
 * @brief Class for implementing an ODEInt::Observer for storing the states at each arc-lengths along the CTR backbone during the numerical integration.
 * @ingroup ctr_utils
 */
class Observer
{
private:
	/**
     * @brief A vector containing all CTR states along the entire backbone.
     */
    std::vector<state_type> &m_states;

    /**
     * @brief A vector containing all discrete arc-length points along the backbone, representing the shape of the CTR.
     */
    std::vector<double> &m_arcLength;

public:
	/**
	 * @brief Removes the default Observer Class constructor
	 */
	Observer() = delete;

	/**
	 * @brief Implements the overloaded constructor for the Observer class.
	 * The Observer is responsible for capturing and storing the values of the state vector as the boost::ODEInt library integrates the set of ODEs modeling the CTR kinematics.
	 *
	 * @param[in,out] states A container that will store the 15-dimensional states along the backbone.
	 * @param[in,out] s A container that will store the discrete arc-length values along the CTR backbone.
	 */
	Observer(std::vector<state_type> &states, std::vector<double> &s) : m_states(states), m_arcLength(s) {}

	/**
	 * @brief Destroys the Observer object.
	 */
	~Observer() = default;

	/**
	 * @brief Functor that overloads the constructor's signature and implements the method for capturing data from Boost::odeInt integrator
	 *
	 * @param[in] states A 15-dimensional static Blaze vector containing the current values of the state vector at the arc-length 's' within the current segment.
	 * @param[in] s The arc-length value corresponding to the current integration step.
	 */
	void operator()(const state_type &states, double s);
};