#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include "Tube.hpp"

/**
 * @brief Class representing Segments between transition points in a Concentric Tube Robot (CTR).
 */
class Segment
{
private:
    /** 
     * @brief Arc-length of each tube transition point.
     */
    std::vector<double> m_S;

    /** 
     * @brief Tubes' bending stiffness -- x, y directions.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_EI;

    /** 
     * @brief Tubes' torsional stiffness -- z direction.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_GJ;

    /** 
     * @brief Tubes' precurvature in the x direction.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_U_x;

    /** 
     * @brief Tubes' precurvature in the y direction.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> m_U_y;

    /** 
     * @brief Arc-length at which precurvature starts in the tubes.
     */
    std::array<double, 3UL> m_len_curv;

    /** 
     * @brief Arc-length at which the tubes terminate (distal-ends).
     */
    std::array<double, 3UL> m_dist_end;

public:
    /**
     * @brief Implements the default constructor for the Segment class.
     */
    Segment();

    /**
     * @brief Implements the overloaded constructor for the Segment class.
     *
     * @param Tb A 3-dimensional std::array containing smart pointers to the three tube objects comprising the CTR assembly.
     * @param beta A 3-dimensional static Blaze vector containing the actuation input values for the linear joints of the CTR.
     */
    Segment(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, const blaze::StaticVector<double, 3UL> &beta);

    /**
     * @brief Implements the copy constructor for the Segment class.
     *
     * @param rhs The source Segment object to copy from.
     */
    Segment(const Segment &rhs);

    /**
     * @brief Implements the move constructor for the Segment class.
     *
     * @param rhs The source Segment object to move from.
     */
    Segment(Segment &&rhs) noexcept;

    /**
     * @brief Destroys the Segment object.
     */
    ~Segment() = default;

    /**
     * @brief Implements the copy assignment operator for the Segment class.
     *
     * @param rhs The source Segment object to copy from.
     * @return A reference to the assigned Segment object.
     */
    Segment &operator=(const Segment &rhs);

    /**
     * @brief Implements the move assignment operator for the Segment class.
     *
     * @param rhs The source Segment object to move from.
     * @return A reference to the assigned Segment object.
     */
    Segment &operator=(Segment &&rhs) noexcept;

    /**
     * @brief Computes the tube transition points and the corresponding parameters at each segment.
     *
     * @param Tb A 3-dimensional std::array containing smart pointers to the three tube objects comprising the CTR assembly.
     * @param beta A 3-dimensional static Blaze vector containing the actuation input values for the linear joints of the CTR.
     */
    void recalculateSegments(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, const blaze::StaticVector<double, 3UL> &beta);

    /**
     * @brief Implements a getter method for retrieving all the transition points currently present in the CTR assembly.
     *
     * @return A std::vector containing the arc-length values (in meters) at which there are tube transition points (where tube ends or where there's a step change in the tube's pre-curvatures).
     */
    std::vector<double> get_S();

    /**
     * @brief Implements a getter method for retrieving the distal ends of all tubes in the CTR assembly.
     *
     * @return A blaze::StaticVector containing the arc-length values (in meters) at which each tube in the CTR assembly terminates.
     */
    blaze::StaticVector<double, 3UL> getDistalEnds();

    /**
     * @brief Implements a getter method for retrieving the bending stiffness of the tubes in all of the tube segments in the CTR assembly.
     *
     * @return A 3xN hybrid Blaze matrix containing the bending stiffness [k_x, k_y, 0] for the tube assembly between the transition points.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_EI();

    /**
     * @brief Implements a getter method for retrieving the torsional stiffness of the tubes in all of the tube segments in the CTR assembly.
     *
     * @return A 3xN hybrid Blaze matrix containing the torsional stiffness [0, 0, k_z] for the tube assembly between the transition points.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_GJ();

    /**
     * @brief Implements a getter method for retrieving the pre-curvature of the tubes along the 'x' direction in all of the tube segments in the CTR assembly.
     *
     * @return A 3xN hybrid Blaze matrix containing the pre-curvatures along the 'x' direction for the tubes in all segments. The 1st, 2nd, and 3rd row of the matrix correspond to the pre-curvatures of the 1st, 2nd, and 3rd tubes, where the 1st tube is the innermost one.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_U_x();

    /**
     * @brief Implements a getter method for retrieving the pre-curvature of the tubes along the 'y' direction in all of the tube segments in the CTR assembly.
     *
     * @return A 3xN hybrid Blaze matrix containing the pre-curvatures along the 'y' direction for the tubes in all segments. The 1st, 2nd, and 3rd row of the matrix correspond to the pre-curvatures of the 1st, 2nd, and 3rd tubes, where the 1st tube is the innermost one.
     */
    blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> get_U_y();

    /**
     * @brief Implements a getter method for retrieving all parameters along all segments in the CTR assembly.
     *
     * @return A Tuple containing: a 3xN Blaze matrix of the bending stiffness (EI), a 3xN Blaze matrix of the torsional stiffness (GJ), a 3xN Blaze matrix of the pre-curvature along the 'x' direction (U_x), a 3xN Blaze matrix of the pre-curvature along the 'y' direction (U_y), and a std::vector with the arc-length (in meters) at which a tube transition occurs (S).
     */
    std::tuple<blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
               blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
               blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
               blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor>,
               std::vector<double>>
    returnParameters();
};