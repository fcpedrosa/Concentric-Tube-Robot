#pragma once

/**
 * @file Tube.hpp
 * @brief A single pre-curved tube (ctr::Tube) and its parameters (ctr::TubeParams).
 */

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
 * All quantities are SI: meters, Pascals, 1/m. A tube consists of a straight
 * transmission section of length `ls` (at the actuator side) followed by a
 * curved section of length `lc` with constant pre-curvature `u_ast`, so its
 * total length is `ls + lc`. See @ref usage_tubes.
 */
struct TubeParams
{
    double OD;                              ///< Outer diameter [m].
    double ID;                              ///< Inner diameter [m].
    double E;                               ///< Young's modulus [Pa].
    double G;                               ///< Shear modulus [Pa].
    double ls;                              ///< Length of the straight transmission section [m].
    double lc;                              ///< Length of the curved section [m].
    blaze::StaticVector<double, 3UL> u_ast; ///< Pre-curvature of the curved section [1/m], tube body frame:
                                            ///< {1/R, 0, 0} for a planar arc of radius R. Only x, y are used.
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
 *
 * @code{.cpp}
 * // Nitinol tube, 1.1 mm OD, straight 120 mm then curved 80 mm with 10 cm radius
 * const ctr::Tube tube{{.OD = 1.10e-3, .ID = 0.97e-3, .E = 65e9, .G = 24.6e9,
 *                       .ls = 120e-3, .lc = 80e-3, .u_ast = {1.0 / 0.1, 0.0, 0.0}}};
 * @endcode
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
     * @throws std::invalid_argument if ID <= 0, OD <= ID, E <= 0, G <= 0,
     *         ls < 0, lc < 0 or ls + lc <= 0 (NaNs are rejected too).
     */
    explicit Tube(const TubeParams &p);

    ~Tube() = default;
    Tube(const Tube &) = default;                ///< Copyable value type.
    Tube(Tube &&) noexcept = default;            ///< Movable.
    Tube &operator=(const Tube &) = default;     ///< Copy-assignable. @return `*this`.
    Tube &operator=(Tube &&) noexcept = default; ///< Move-assignable. @return `*this`.

    // ─── Getters ─────────────────────────────────────────────────────────────

    /**
     * @brief Returns all physical parameters packed in a named aggregate.
     * @return A copy of the parameters the tube was built from.
     */
    [[nodiscard]] TubeParams parameters() const noexcept;

    /**
     * @brief Returns the total tube length.
     * @return ls + lc [m].
     */
    [[nodiscard]] double getTubeLength() const noexcept;

    /**
     * @brief Returns the straight-section length.
     * @return ls [m].
     */
    [[nodiscard]] double getStraightLen() const noexcept;

    /**
     * @brief Returns the curved-section length.
     * @return lc [m].
     */
    [[nodiscard]] double getCurvLen() const noexcept;

    /**
     * @brief Returns the pre-curvature vector.
     * @return u* [1/m], tube body frame.
     */
    [[nodiscard]] blaze::StaticVector<double, 3UL> get_u_ast() const noexcept;

    /**
     * @brief Returns one component of the pre-curvature vector.
     * @pre `id < 3` (not checked).
     * @param id 0-based component index: 0 = x, 1 = y, 2 = z.
     * @return u*[id] [1/m].
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
