#pragma once

#include <blaze/Math.h>
#include <vector>

typedef blaze::StaticVector<double, 15UL> state_type;

/**
 * @brief Class for implementing an ODEInt::Observer for storing the states at each arc-lengths along the CTR backbone during the numerical integration.
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
	 * @param states A 15-dimensional static Blaze vector containing the current values of the state vector at the arc-length 's' within the current segment.
	 * @param s A std::vector containing the set discrete arc-length values along along the CTR backbone.
	 */
	Observer(std::vector<state_type> &states, std::vector<double> &s) : m_states(states), m_arcLength(s) {}

	/**
	 * @brief Destroys the Observer object.
	 */
	~Observer() = default;

	/**
	 * @brief Functor that overloads the constructor's signature and implements the method for capturing data from Boost::odeInt integrator
	 *
	 * @param states A 15-dimensional static Blaze vector containing the current values of the state vector at the arc-length 's' within the current segment.
	 * @param s A std::vector containing the set discrete arc-length values along along the CTR backbone.
	 */
	void operator()(const state_type &states, double s);
};