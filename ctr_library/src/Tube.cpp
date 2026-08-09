#include "Tube.hpp"
#include <cassert>

namespace ctr {

// ─── Constructors ─────────────────────────────────────────────────────────────

Tube::Tube()
    : m_OD(0.0), m_ID(0.0), m_E(0.0), m_G(0.0),
      m_ls(0.0), m_lc(0.0), m_u_ast(0.0)
{}

Tube::Tube(double OD, double ID, double E, double G, double ls, double lc,
           const blaze::StaticVector<double, 3UL> &u_ast)
    : m_OD(OD), m_ID(ID), m_E(E), m_G(G), m_ls(ls), m_lc(lc), m_u_ast(u_ast)
{
    assert(OD > ID   && "Tube: outer diameter must exceed inner diameter");
    assert(ID > 0.0  && "Tube: inner diameter must be positive");
    assert(E  > 0.0  && "Tube: Young's modulus must be positive");
    assert(G  > 0.0  && "Tube: shear modulus must be positive");
    assert((ls + lc) > 0.0 && "Tube: total tube length must be positive");
}

// ─── Getters ─────────────────────────────────────────────────────────────────

TubeParams Tube::getParameters() const noexcept
{
    return {m_OD, m_ID, m_E, m_G, m_ls, m_lc, m_u_ast};
}

double Tube::getTubeLength() const noexcept
{
    return m_ls + m_lc;
}

double Tube::getStraightLen() const noexcept
{
    return m_ls;
}

double Tube::getCurvLen() const noexcept
{
    return m_lc;
}

blaze::StaticVector<double, 3UL> Tube::get_u_ast() const noexcept
{
    return m_u_ast;
}

double Tube::get_u_ast(std::size_t id) const noexcept
{
    assert(id <= 1UL && "get_u_ast: id must be 0 (x) or 1 (y)");
    return m_u_ast[id];
}

double Tube::getK(Stiffness s) const noexcept
{
    const double csf = crossSectionFactor();
    switch (s)
    {
    case Stiffness::Bending: return m_E * pi_64 * csf;
    case Stiffness::Torsion: return m_G * pi_32 * csf;
    }
    return 0.0; // unreachable — satisfies compiler
}

blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>>
Tube::getK_Matrix() const noexcept
{
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K;
    const double EI = getK(Stiffness::Bending);
    const double GJ = getK(Stiffness::Torsion);
    K(0UL, 0UL) = K(1UL, 1UL) = EI;
    K(2UL, 2UL)                = GJ;
    return K;
}

// ─── Setters ─────────────────────────────────────────────────────────────────

void Tube::setYoungModulus(double E) noexcept
{
    m_E = E;
    // Stiffness is computed on demand — no cache to synchronise.
}

void Tube::setShearModulus(double G) noexcept
{
    m_G = G;
}

void Tube::setK(double EI, double GJ) noexcept
{
    // Back-compute E and G from the given stiffness values.
    const double csf = crossSectionFactor();
    if (const double I = pi_64 * csf; I > 0.0) m_E = EI / I;
    if (const double J = pi_32 * csf; J > 0.0) m_G = GJ / J;
}

void Tube::setBendingK(double EI) noexcept
{
    if (const double I = pi_64 * crossSectionFactor(); I > 0.0)
        m_E = EI / I;
}

void Tube::set_u_ast(const blaze::StaticVector<double, 3UL> &u_ast)
{
    m_u_ast = u_ast;
}

void Tube::set_u_ast(std::size_t id, double u)
{
    assert(id <= 1UL && "set_u_ast: id must be 0 (x) or 1 (y)");
    m_u_ast[id] = u;
}

void Tube::setStraightLen(double ls) noexcept
{
    m_ls = (ls > 0.0) ? ls : 0.0;
}

void Tube::setCurvLen(double lc) noexcept
{
    m_lc = (lc > 0.0) ? lc : 0.0;
}

} // namespace ctr
