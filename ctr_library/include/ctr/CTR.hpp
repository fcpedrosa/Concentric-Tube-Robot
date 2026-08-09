#pragma once

#include "ctr/Types.hpp"
#include "ctr/Tube.hpp"
#include "ctr/Segment.hpp"
#include "ctr/ODESystem.hpp"

#include <array>
#include <memory>
#include <span>
#include <vector>

namespace ctr
{

class BVPSolver;       // private strategy hierarchy (implementation detail, defined in src/)
class ShootingProblem; // facade granting solvers access to the FK internals

/**
 * @brief A three-tube Concentric Tube Robot (CTR).
 *
 * Responsibilities:
 *  - Maintain robot state (tube assembly, joint configuration, external loads).
 *  - Perform forward kinematics by solving the Cosserat-model boundary value
 *    problem with a shooting method (Rucker, Jones & Webster, IEEE T-RO 2010).
 *  - Solve inverse kinematics (tip position control).
 *
 * Conventions (SI units throughout):
 *  - Joint configuration q = [β₁, β₂, β₃, α₁, α₂, α₃]: linear retractions β
 *    in meters (non-positive), axial rotations α in radians.
 *  - Tubes are ordered innermost-first everywhere.
 *  - All returned positions/shapes are in meters, in the global frame at the
 *    robot base (s = 0).
 */
class CTR
{
  public:
    CTR() = delete;

    /**
     * @brief Constructs a CTR robot.
     *
     * @param tubes        The three component tubes, innermost first.
     * @param q            Initial 6-DOF actuation vector [β₁..β₃ (m), α₁..α₃ (rad)].
     * @param bvpTolerance Convergence tolerance for the BVP residue (L∞ norm).
     * @param method       Root-finding method for the shooting BVP.
     */
    CTR(std::array<Tube, NUM_TUBES> tubes, const blaze::StaticVector<double, 6UL> &q, double bvpTolerance,
        RootFindingMethod method = RootFindingMethod::ModifiedNewtonRaphson);

    CTR(const CTR &rhs);
    CTR &operator=(const CTR &rhs);
    CTR(CTR &&rhs) noexcept;
    CTR &operator=(CTR &&rhs) noexcept;
    ~CTR();

    // ─── Kinematics ──────────────────────────────────────────────────────────

    /**
     * @brief Actuates the CTR to a joint configuration and solves the corresponding BVP.
     *
     * @param q         Target 6-DOF joint configuration [β (m), α (rad)].
     * @param initGuess Proximal boundary-condition guess (updated in place on
     *                  convergence — pass the previous solution to warm-start).
     * @return FKResult describing convergence, iterations, and final residue.
     */
    [[nodiscard]] FKResult actuate(const blaze::StaticVector<double, 6UL> &q, bvp_type &initGuess);

    /**
     * @brief Inverse kinematics: steers the tip to a Cartesian target position.
     *
     * @param target    Desired 3D tip position [m], global frame.
     * @param posTol    Position tolerance [m]; success means ||tip − target|| ≤ posTol.
     * @param initGuess Proximal boundary-condition guess (updated in place; warm-started across calls).
     * @return IKResult; IKResult::converged reports target attainment (not merely BVP convergence).
     */
    [[nodiscard]] IKResult solveIK(const blaze::StaticVector<double, 3UL> &target, double posTol,
                                   bvp_type &initGuess);

    /**
     * @brief Finite-difference kinematic Jacobian ∂tip/∂q (3×6).
     *
     * Exposed for research use; solveIK computes what it needs internally.
     *
     * @param initGuess Current proximal boundary-condition guess.
     * @param tipPos    Current tip position [m].
     * @return 3×6 kinematic Jacobian.
     */
    [[nodiscard]] Mat<3UL, 6UL> kinematicJacobian(const bvp_type &initGuess,
                                                  const blaze::StaticVector<double, 3UL> &tipPos);

    // ─── Shape access (all meters) ───────────────────────────────────────────

    /// A polyline of 3D backbone points [m].
    using Points = std::vector<blaze::StaticVector<double, 3UL>>;

    /**
     * @brief Returns the centerline of each tube as a polyline of 3D points [m].
     * @return Array indexed innermost-first; tube 1 spans the full backbone.
     */
    [[nodiscard]] std::array<Points, NUM_TUBES> tubeShapes() const;

    /** @brief Returns the full backbone centerline as a polyline of 3D points [m]. */
    [[nodiscard]] Points backboneShape() const;

    /** @brief Returns the tip (distal end) position [m]. */
    [[nodiscard]] blaze::StaticVector<double, 3UL> tipPosition() const;

    /** @brief Returns the arc-length of each tube's distal end [m]. */
    [[nodiscard]] blaze::StaticVector<double, NUM_TUBES> distalEnds() const;

    // ─── Observers ───────────────────────────────────────────────────────────

    /** @brief Returns the tube assembly (innermost first). */
    [[nodiscard]] const std::array<Tube, NUM_TUBES> &tubes() const noexcept;

    /** @brief Returns the current joint configuration [β (m), α (rad)]. */
    [[nodiscard]] blaze::StaticVector<double, 6UL> configuration() const noexcept;

    /** @brief Returns the linear actuation values [β₁, β₂, β₃] (m). */
    [[nodiscard]] blaze::StaticVector<double, NUM_TUBES> beta() const noexcept;

    /** @brief Returns the BVP convergence tolerance. */
    [[nodiscard]] double tolerance() const noexcept { return m_accuracy; }

    /** @brief Returns the fixed arc-length integration step [m]. */
    [[nodiscard]] double integrationStep() const noexcept { return m_ds; }

    /** @brief Full ODE state at each recorded arc-length sample (advanced). */
    [[nodiscard]] std::span<const state_type> states() const noexcept;

    /** @brief Arc-length values corresponding to states() (advanced). */
    [[nodiscard]] std::span<const double> arcLengthSamples() const noexcept;

    // ─── Modifiers ───────────────────────────────────────────────────────────

    /** @brief Sets the joint configuration without actuating the CTR. */
    void setConfiguration(const blaze::StaticVector<double, 6UL> &q);

    /** @brief Switches the BVP root-finding strategy. */
    void setBVPMethod(RootFindingMethod method);

    /** @brief Sets the external point moment at the CTR's distal end [N·m], global frame. */
    void setDistalMoment(const blaze::StaticVector<double, 3UL> &moment);

    /** @brief Sets the external point force at the CTR's distal end [N], global frame. */
    void setDistalForce(const blaze::StaticVector<double, 3UL> &force);

    /** @brief Replaces one tube (innermost-first index) and recomputes the segmentation. */
    void setTube(std::size_t idx, Tube tube);

    /**
     * @brief Sets the fixed arc-length integration step [m] (default 1 mm).
     *
     * Integration is deliberately fixed-step and deterministic (the
     * finite-difference Jacobians depend on it); this knob trades accuracy
     * for speed uniformly. Values in [1e-5, 1e-2] are accepted.
     */
    void setIntegrationStep(double ds);

  private:
    friend class ShootingProblem; // the only external access path to the FK internals

    // ─── FK internals ─────────────────────────────────────────────────────────

    /** Resets integration buffers and initial conditions prior to an ODE shot. */
    void reset(const bvp_type &initGuess);

    /** Integrates the backbone ODE once; returns the distal boundary-condition residue. */
    [[nodiscard]] bvp_type residual(const bvp_type &initGuess);

    /** Finite-difference Jacobian of the BVP residue w.r.t. the shooting variables (5×5). */
    [[nodiscard]] Mat<BVP_DIM, BVP_DIM> bvpJacobian(const bvp_type &initGuess, const bvp_type &residue);

    /** Lazily recreates the solver if this object was moved from. */
    void ensureSolver();

    /** Returns the linear actuation sub-vector [β₁, β₂, β₃] as a view into m_q. */
    [[nodiscard]] decltype(auto) betaView() const noexcept { return blaze::subvector<0UL, NUM_TUBES>(m_q); }

    // ─── Physical model ───────────────────────────────────────────────────────

    /** The tube assembly (innermost first). */
    std::array<Tube, NUM_TUBES> m_tubes;

    /** Joint actuation values [β₁, β₂, β₃, α₁, α₂, α₃]. */
    blaze::StaticVector<double, 6UL> m_q;

    /** Initial twist angles for all tubes at s = 0 (scratch, set by reset()). */
    blaze::StaticVector<double, NUM_TUBES> m_theta_0;

    /** Initial orientation of the local frame at s = 0 (quaternion; scratch, set by reset()). */
    blaze::StaticVector<double, 4UL> m_h_0;

    /** External point force at the distal end [N]. */
    blaze::StaticVector<double, 3UL> m_wf;

    /** External point moment at the distal end [N·m]. */
    blaze::StaticVector<double, 3UL> m_wm;

    // ─── Computation components ───────────────────────────────────────────────

    /** BVP convergence tolerance. */
    double m_accuracy;

    /** Fixed arc-length integration step [m]. */
    double m_ds = 1.0e-3;

    /** Active root-finding method — stored so copy/move preserve the solver type. */
    RootFindingMethod m_method;

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
};

} // namespace ctr
