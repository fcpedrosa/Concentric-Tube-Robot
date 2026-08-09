// ***************************************************************************************** //
// *  This file is part of NIH_CSTAR, a CTR kinematics library for Concentric Tube Robots  * //
// *
// * //
// *  ----------- # Copyright (C) 2021 Filipe C. Pedrosa <fpedrosa@uwo.ca> # -----------   * //
// *
// * //
// *  Project developed under the supervision of Prof Dr Rajni Patel <rvpatel@uwo.ca>	   * //
// *			  CSTAR (Canadian Surgical Technologies & Advanced Robotics)			   * //
// *				   Western University, London, ON, Canada
// * //
// ***************************************************************************************** //

#pragma once

#include "Tube.hpp"
#include "Segment.hpp"
#include "ODESystem.hpp"
#include "BVPSolver.hpp"
#include "mathOperations.hpp"
#include "boostBlazeAlgebra.hpp"
#include <array>
#include <memory>
#include <span>
#include <tuple>
#include <vector>

namespace ctr
{

/**
 * @brief Class implementing a three-tube Concentric Tube Robot (CTR).
 *
 * Responsibilities:
 *  - Maintain robot state (tube assembly, joint configuration).
 *  - Perform forward kinematics by integrating the backbone ODE system.
 *  - Compute the BVP and kinematic Jacobians.
 *  - Delegate BVP root-finding to a pluggable BVPSolver strategy.
 *  - Solve inverse kinematics (position control) via resolved-rates motion.
 */
class CTR
{
  public:
    CTR() = delete;

    /**
     * @brief Constructs a CTR object.
     *
     * @param Tb     Array of shared pointers to the three Tube objects (innermost first).
     * @param q      6-DOF actuation vector [beta_1..3, alpha_1..3].
     * @param Tol    Convergence tolerance for the BVP root-finding.
     * @param method Root-finding method to use.
     */
    CTR(const std::array<std::shared_ptr<Tube>, NUM_TUBES> &Tb, const blaze::StaticVector<double, 6UL> &q, double Tol,
        mathOp::rootFindingMethod method);

    CTR(const CTR &rhs);
    CTR(CTR &&rhs) noexcept = default;
    ~CTR() = default;
    CTR &operator=(const CTR &rhs);
    CTR &operator=(CTR &&rhs) noexcept = default;

    // ─── Forward kinematics ──────────────────────────────────────────────────

    /**
     * @brief Resets the integration buffers and initial conditions prior to ODE solving.
     * @param initGuess Proximal boundary condition guess [mb_x(0), mb_y(0), u1z(0), u2z(0), u3z(0)].
     */
    void reset(const bvp_type &initGuess);

    /**
     * @brief Integrates the CTR ODE system and returns the BVP residue.
     * @param initGuess Proximal boundary condition guess.
     * @return 5-dimensional residue vector (violation of distal boundary conditions).
     */
    [[nodiscard]] bvp_type ODESolver(const bvp_type &initGuess);

    /**
     * @brief Computes the finite-difference Jacobian of the BVP residue w.r.t. the initial guess.
     * @param initGuess Current proximal boundary condition guess.
     * @param residue   BVP residue at initGuess.
     * @return 5×5 BVP Jacobian matrix.
     */
    [[nodiscard]] blaze::StaticMatrix<double, BVP_DIM, BVP_DIM> jac_BVP(const bvp_type &initGuess,
                                                                        const bvp_type &residue);

    /**
     * @brief Computes the finite-difference kinematic Jacobian (tip velocity w.r.t. joint velocities).
     * @param initGuess Current proximal boundary condition guess.
     * @param tipPos    Current distal-end (tip) position.
     * @return 3×6 kinematic Jacobian matrix.
     */
    [[nodiscard]] blaze::StaticMatrix<double, 3UL, 6UL> jacobian(const bvp_type &initGuess,
                                                                 const blaze::StaticVector<double, 3UL> &tipPos);

    // ─── Actuation & inverse kinematics ─────────────────────────────────────

    /**
     * @brief Actuates the CTR to a given joint configuration and solves the corresponding BVP.
     * @param initGuess Proximal boundary condition guess (updated in place on convergence).
     * @param q_input   Target 6-DOF joint configuration.
     * @return true if the BVP converged, false otherwise.
     */
    [[nodiscard]] bool actuate_CTR(bvp_type &initGuess, const blaze::StaticVector<double, 6UL> &q_input);

    /**
     * @brief Resolved-rates position controller. Steers the tip to a target position.
     * @param initGuess Proximal boundary condition guess.
     * @param target    Desired 3D tip position.
     * @param posTol    Convergence tolerance in meters.
     * @return true if the target was reached within posTol, false otherwise.
     */
    [[nodiscard]] bool posCTRL(bvp_type &initGuess, const blaze::StaticVector<double, 3UL> &target, double posTol);

    // ─── Getters ─────────────────────────────────────────────────────────────

    [[nodiscard]] std::array<std::shared_ptr<Tube>, NUM_TUBES> getTubes() const;
    [[nodiscard]] blaze::StaticVector<double, 3UL> getBeta() const;
    [[nodiscard]] blaze::StaticVector<double, 6UL> getConfiguration() const;
    [[nodiscard]] blaze::StaticVector<double, 3UL> getTipPos() const;
    [[nodiscard]] blaze::StaticVector<double, NUM_TUBES> getDistalEnds() const;

    /** @return Tuple of three 3×N matrices giving x,y,z coordinates of each tube's centerline. */
    [[nodiscard]] std::tuple<blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>,
                             blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>,
                             blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>>
    getTubeShapes() const;

    /** @return Tuple of three std::vectors with x,y,z coordinates of the innermost tube. */
    [[nodiscard]] std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> getShape() const;

    [[nodiscard]] std::span<const state_type> states() const noexcept;
    [[nodiscard]] std::span<const double> arcLengthSamples() const noexcept;

    /** @return The BVP convergence tolerance. */
    [[nodiscard]] double accuracy() const noexcept { return m_accuracy; }

    // ─── Setters ─────────────────────────────────────────────────────────────

    /** @brief Sets the joint configuration without actuating the CTR. */
    void setConfiguration(const blaze::StaticVector<double, 6UL> &q);

    /** @brief Switches the BVP root-finding strategy. */
    void setBVPMethod(mathOp::rootFindingMethod mthd);

    /** @brief Sets the external point moment at the CTR's distal end. */
    void setDistalMoment(const blaze::StaticVector<double, 3UL> &moment);

    /** @brief Sets the external point force at the CTR's distal end. */
    void setDistalForce(const blaze::StaticVector<double, 3UL> &force);

  private:
    // ─── Physical model ───────────────────────────────────────────────────────

    /** Array of shared pointers to the Tube objects comprising the CTR. */
    std::array<std::shared_ptr<Tube>, NUM_TUBES> m_Tubes;

    /** Joint actuation values [beta_1, beta_2, beta_3, alpha_1, alpha_2, alpha_3]. */
    blaze::StaticVector<double, 6UL> m_q;

    /** Initial twist angles for all tubes at s = 0. */
    blaze::StaticVector<double, NUM_TUBES> m_theta_0;

    /** Initial orientation of the local frame at s = 0 (quaternion). */
    blaze::StaticVector<double, 4UL> m_h_0;

    /** External point force at the distal end. */
    blaze::StaticVector<double, 3UL> m_wf;

    /** External point moment at the distal end. */
    blaze::StaticVector<double, 3UL> m_wm;

    // ─── Computation components ───────────────────────────────────────────────

    /** BVP convergence tolerance. */
    double m_accuracy;

    /** Active root-finding method — stored so copy/move preserve the solver type. */
    mathOp::rootFindingMethod m_method;

    /** Tube transition points and per-segment kinematic parameters. */
    Segment m_segment;

    /** ODE system implementing the CTR backbone kinematics. */
    ODESystem m_stateEquations;

    /** ODE integration state history (one entry per arc-length step). */
    std::vector<state_type> m_y;

    /** Arc-length values corresponding to each entry in m_y. */
    std::vector<double> m_s;

    /** Pluggable BVP root-finding strategy. */
    std::unique_ptr<BVPSolver> m_solver;

    // ─── Helpers ─────────────────────────────────────────────────────────────

    /** Returns the linear actuation sub-vector [beta_1, beta_2, beta_3] from m_q. */
    [[nodiscard]] decltype(auto) beta() const noexcept { return blaze::subvector<0UL, NUM_TUBES>(m_q); }

    /**
     * @brief Extracts raw read-only pointers to the three Tube objects.
     *
     * Used to pass to Segment without transferring ownership or requiring
     * shared_ptr in the Segment interface.
     */
    [[nodiscard]] std::array<const Tube *, NUM_TUBES> rawTubes() const noexcept
    {
        return {m_Tubes[0].get(), m_Tubes[1].get(), m_Tubes[2].get()};
    }
};

} // namespace ctr
