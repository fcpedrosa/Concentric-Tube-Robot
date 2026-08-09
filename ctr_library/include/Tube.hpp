#pragma once

#include <cmath>
#include <numbers>
#include <cassert>
#include <blaze/Math.h>

namespace ctr {

// ─── Tube stiffness selector ───────────────────────────────────────────────────

/**
 * @brief Selects which stiffness component to retrieve from a Tube.
 *
 * Replaces the error-prone integer index previously passed to getK().
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
 * Returned by Tube::getParameters(). A named struct eliminates the
 * positional confusion of the old std::tuple return value.
 */
struct TubeParams
{
    double OD;   ///< Outer diameter [m].
    double ID;   ///< Inner diameter [m].
    double E;    ///< Young's modulus [Pa].
    double G;    ///< Shear modulus [Pa].
    double ls;   ///< Length of the straight transmission section [m].
    double lc;   ///< Length of the curved section [m].
    blaze::StaticVector<double, 3UL> u_ast; ///< Pre-curvature vector [1/m].
};

// ─── Tube ─────────────────────────────────────────────────────────────────────

/**
 * @brief Represents a single tube in the concentric arrangement of a CTR robot.
 *
 * Stores raw geometry (OD, ID) and material constants (E, G).
 * Derived quantities (I, J, EI, GJ) are computed on demand by getK() to
 * eliminate the synchronisation hazard that arose from caching them alongside
 * the raw parameters.
 */
class Tube
{
private:
    static constexpr double pi_64 = std::numbers::pi / 64.0; ///< π/64 — factor in I.
    static constexpr double pi_32 = std::numbers::pi / 32.0; ///< π/32 — factor in J.

    double m_OD;  ///< Outer diameter [m].
    double m_ID;  ///< Inner diameter [m].
    double m_E;   ///< Young's modulus [Pa].
    double m_G;   ///< Shear modulus [Pa].
    double m_ls;  ///< Straight-section length [m].
    double m_lc;  ///< Curved-section length [m].
    blaze::StaticVector<double, 3UL> m_u_ast; ///< Pre-curvature vector [1/m].

    /// Returns OD^4 − ID^4, shared by bending and torsional stiffness formulas.
    [[nodiscard]] double crossSectionFactor() const noexcept
    {
        const double od2 = m_OD * m_OD;
        const double id2 = m_ID * m_ID;
        return od2 * od2 - id2 * id2;
    }

public:
    // ─── Constructors / special members ──────────────────────────────────────

    /**
     * @brief Default constructor — creates an invalid (all-zero) tube.
     * @note Call isValid() before using any computed property.
     */
    Tube();

    /**
     * @brief Constructs a fully specified Tube.
     *
     * @pre OD > ID > 0, E > 0, G > 0, (ls + lc) > 0.
     *
     * @param OD    Outer diameter [m].
     * @param ID    Inner diameter [m].
     * @param E     Young's modulus [Pa].
     * @param G     Shear modulus [Pa].
     * @param ls    Straight-section length [m].
     * @param lc    Curved-section length [m].
     * @param u_ast Pre-curvature vector [1/m].
     */
    Tube(double OD, double ID, double E, double G, double ls, double lc,
         const blaze::StaticVector<double, 3UL> &u_ast);

    ~Tube()                           = default;
    Tube(const Tube &)                = default;
    Tube(Tube &&) noexcept            = default;
    Tube &operator=(const Tube &)     = default;
    Tube &operator=(Tube &&) noexcept = default;

    // ─── Validation ──────────────────────────────────────────────────────────

    /**
     * @brief Returns true if the tube has been initialised with valid geometry.
     */
    [[nodiscard]] bool isValid() const noexcept
    {
        return m_OD > 0.0 && m_ID > 0.0 && m_OD > m_ID && (m_ls + m_lc) > 0.0;
    }

    // ─── Getters ─────────────────────────────────────────────────────────────

    /**
     * @brief Returns all physical parameters packed in a named aggregate.
     */
    [[nodiscard]] TubeParams getParameters() const noexcept;

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
     */
    [[nodiscard]] double getK(Stiffness s) const noexcept;

    /**
     * @brief Returns the full 3×3 diagonal stiffness matrix [EI, EI, GJ], on demand.
     */
    [[nodiscard]] blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>>
    getK_Matrix() const noexcept;

    // ─── Setters ─────────────────────────────────────────────────────────────

    /**
     * @brief Updates Young's modulus [Pa]. getK(Bending) reflects the new value immediately.
     */
    void setYoungModulus(double E) noexcept;

    /**
     * @brief Updates the shear modulus [Pa]. getK(Torsion) reflects the new value immediately.
     */
    void setShearModulus(double G) noexcept;

    /**
     * @brief Sets both stiffness values by back-computing E and G from the cross-section geometry.
     * @param EI New bending stiffness [N·m²].
     * @param GJ New torsional stiffness [N·m²].
     */
    void setK(double EI, double GJ) noexcept;

    /**
     * @brief Sets only the bending stiffness by back-computing E.
     * @param EI New bending stiffness [N·m²].
     */
    void setBendingK(double EI) noexcept;

    /** @brief Updates the pre-curvature vector [1/m]. */
    void set_u_ast(const blaze::StaticVector<double, 3UL> &u_ast);

    /**
     * @brief Updates the pre-curvature along one direction.
     * @param id 0-based: 0 = x, 1 = y.
     * @param u  New scalar pre-curvature [1/m].
     */
    void set_u_ast(std::size_t id, double u);

    /** @brief Updates the straight-section length [m]. */
    void setStraightLen(double ls) noexcept;

    /** @brief Updates the curved-section length [m]. */
    void setCurvLen(double lc) noexcept;
};

} // namespace ctr
