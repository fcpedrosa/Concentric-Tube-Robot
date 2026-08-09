#pragma once

#include <vector>
#include <array>
#include "ctr/Tube.hpp"
#include "ctr/Types.hpp"

namespace ctr
{

class CTR; // forward declaration — CTR is the only class that may recalculate segments

// ─── Segment ───────────────────────────────────────────────────────────────────

/**
 * @brief Computes and stores the per-segment kinematic parameters of a CTR.
 *
 * A CTR backbone is divided into segments by the arc-length positions at which
 * tubes start or end. This class finds those transition points and caches the
 * per-segment stiffness and pre-curvature matrices needed by the ODE integrator.
 *
 * recalculateSegments() is private: only the owning CTR object is permitted to
 * update the segment decomposition (via the friend declaration).
 */
class Segment
{
    friend class CTR;

  private:
    std::vector<double> m_S; ///< Arc-length transition points.

    blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> m_EI;  ///< Bending stiffness (3 × N).
    blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> m_GJ;  ///< Torsional stiffness (3 × N).
    blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> m_U_x; ///< Pre-curvature x (3 × N).
    blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> m_U_y; ///< Pre-curvature y (3 × N).

    blaze::StaticVector<double, NUM_TUBES> m_len_curv; ///< Arc-length at which each tube's curved section starts.
    blaze::StaticVector<double, NUM_TUBES> m_dist_end; ///< Arc-length at each tube's distal end.

    /**
     * @brief Recomputes all tube transition points and per-segment parameters.
     *
     * Must be called (by CTR) whenever the linear actuation inputs change.
     *
     * @param tubes The NUM_TUBES Tube objects (innermost first).
     * @param beta  Linear actuation inputs [beta_1, beta_2, beta_3].
     */
    void recalculateSegments(const std::array<Tube, NUM_TUBES> &tubes,
                             const blaze::StaticVector<double, NUM_TUBES> &beta);

  public:
    Segment() = default;

    /**
     * @brief Constructs and immediately computes the segment decomposition.
     *
     * @param tubes The NUM_TUBES Tube objects (innermost first).
     * @param beta  Linear actuation inputs [beta_1, beta_2, beta_3].
     */
    Segment(const std::array<Tube, NUM_TUBES> &tubes, const blaze::StaticVector<double, NUM_TUBES> &beta);

    Segment(const Segment &) = default;                ///< Copyable.
    Segment(Segment &&) noexcept = default;             ///< Movable.
    ~Segment() = default;
    Segment &operator=(const Segment &) = default;      ///< Copy-assignable.
    Segment &operator=(Segment &&) noexcept = default;  ///< Move-assignable.

    // ─── Getters (cheap, by const reference) ─────────────────────────────────

    /** @brief Returns the arc-length transition points along the CTR backbone. */
    [[nodiscard]] const std::vector<double> &getTransitionPoints() const noexcept;

    /** @brief Returns the arc-lengths at each tube's distal end. */
    [[nodiscard]] const blaze::StaticVector<double, NUM_TUBES> &getDistalEnds() const noexcept;

    /** @brief Returns the bending stiffness matrix over all segments (3 × N). */
    [[nodiscard]] const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &
    get_EI() const noexcept;

    /** @brief Returns the torsional stiffness matrix over all segments (3 × N). */
    [[nodiscard]] const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &
    get_GJ() const noexcept;

    /** @brief Returns the x-direction pre-curvature matrix over all segments (3 × N). */
    [[nodiscard]] const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &
    get_U_x() const noexcept;

    /** @brief Returns the y-direction pre-curvature matrix over all segments (3 × N). */
    [[nodiscard]] const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &
    get_U_y() const noexcept;
};

} // namespace ctr
