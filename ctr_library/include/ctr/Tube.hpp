#pragma once

#include <cmath>
#include <numbers>
#include <blaze/Math.h>

namespace ctr
{

// ─── Tube stiffness selector ───────────────────────────────────────────────────

/**
 * @brief Selects which stiffness component to retrieve from a Tube.
 */
enum class Stiffness
{
    Bending, ///< Bending stiffness EI (identical along x and y).
    Torsion  ///< Torsional stiffness GJ (along z).
};

// ─── Plain aggregate for tube parameters ──────────────────────────────────────

/**
 * @brief Aggregates all physical parameters of a Tube.
 *
 * Used both to construct a Tube (conveniently with designated initializers)
 * and as the return value of Tube::parameters().
 *
 * All quantities are SI: meters, Pascals, 1/m.
 */
struct TubeParams
{
    double OD;                              ///< Outer diameter [m].
    double ID;                              ///< Inner diameter [m].
    double E;                               ///< Young's modulus [Pa].
    double G;                               ///< Shear modulus [Pa].
    double ls;                              ///< Length of the straight transmission section [m].
    double lc;                              ///< Length of the curved section [m].
    blaze::StaticVector<double, 3UL> u_ast; ///< Pre-curvature vector [1/m] (z component must be 0).
};

// ─── Tube ─────────────────────────────────────────────────────────────────────

/**
 * @brief Represents a single tube in the concentric arrangement of a CTR robot.
 *
 * An immutable value type: geometry and material are validated once at
 * construction (throwing std::invalid_argument on bad input), after which a
 * Tube that exists is always valid. To "modify" a tube, edit a TubeParams
 * copy from parameters() and construct a new Tube.
 *
 * Derived quantities (I, J, EI, GJ) are computed on demand by getK().
 */
class Tube
{
  private:
    static constexpr double pi_64 = std::numbers::pi / 64.0; ///< π/64 — factor in I.
    static constexpr double pi_32 = std::numbers::pi / 32.0; ///< π/32 — factor in J.

    double m_OD;                              ///< Outer diameter [m].
    double m_ID;                              ///< Inner diameter [m].
    double m_E;                               ///< Young's modulus [Pa].
    double m_G;                               ///< Shear modulus [Pa].
    double m_ls;                              ///< Straight-section length [m].
    double m_lc;                              ///< Curved-section length [m].
    blaze::StaticVector<double, 3UL> m_u_ast; ///< Pre-curvature vector [1/m].

    /// Returns OD^4 − ID^4, shared by bending and torsional stiffness formulas.
    [[nodiscard]] double crossSectionFactor() const noexcept
    {
        const double od2 = m_OD * m_OD;
        const double id2 = m_ID * m_ID;
        return od2 * od2 - id2 * id2;
    }

  public:
    Tube() = delete;

    /**
     * @brief Constructs a fully specified Tube.
     *
     * @param p Physical parameters (SI units).
     * @throws std::invalid_argument if OD <= ID, ID <= 0, E <= 0, G <= 0,
     *         or (ls + lc) <= 0.
     */
    explicit Tube(const TubeParams &p);

    ~Tube() = default;
    Tube(const Tube &) = default;
    Tube(Tube &&) noexcept = default;
    Tube &operator=(const Tube &) = default;
    Tube &operator=(Tube &&) noexcept = default;

    // ─── Getters ─────────────────────────────────────────────────────────────

    /** @brief Returns all physical parameters packed in a named aggregate. */
    [[nodiscard]] TubeParams parameters() const noexcept;

    /** @brief Returns the total tube length (straight + curved) [m]. */
    [[nodiscard]] double getTubeLength() const noexcept;

    /** @brief Returns the straight-section length [m]. */
    [[nodiscard]] double getStraightLen() const noexcept;

    /** @brief Returns the curved-section length [m]. */
    [[nodiscard]] double getCurvLen() const noexcept;

    /** @brief Returns the pre-curvature vector [1/m]. */
    [[nodiscard]] blaze::StaticVector<double, 3UL> get_u_ast() const noexcept;

    /**
     * @brief Returns the pre-curvature along one direction.
     * @param id 0-based: 0 = x, 1 = y.
     */
    [[nodiscard]] double get_u_ast(std::size_t id) const noexcept;

    /**
     * @brief Returns the requested stiffness coefficient, computed on demand.
     *
     * - Bending: EI = E × (π/64) × (OD⁴ − ID⁴)
     * - Torsion: GJ = G × (π/32) × (OD⁴ − ID⁴)
     *
     * @param s Stiffness::Bending or Stiffness::Torsion.
     * @return The stiffness [N·m²].
     */
    [[nodiscard]] double getK(Stiffness s) const noexcept;
};

} // namespace ctr
