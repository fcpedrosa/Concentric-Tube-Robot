#include "ctr/Tube.hpp"
#include <stdexcept>

namespace ctr
{

// ─── Constructor ──────────────────────────────────────────────────────────────

Tube::Tube(const TubeParams &p)
    : m_OD(p.OD), m_ID(p.ID), m_E(p.E), m_G(p.G), m_ls(p.ls), m_lc(p.lc), m_u_ast(p.u_ast)
{
    // Construction is a cold path — validate unconditionally (also in Release).
    if (!(p.ID > 0.0))
        throw std::invalid_argument("Tube: inner diameter must be positive");
    if (!(p.OD > p.ID))
        throw std::invalid_argument("Tube: outer diameter must exceed inner diameter");
    if (!(p.E > 0.0))
        throw std::invalid_argument("Tube: Young's modulus must be positive");
    if (!(p.G > 0.0))
        throw std::invalid_argument("Tube: shear modulus must be positive");
    if (!(p.ls >= 0.0) || !(p.lc >= 0.0) || !(p.ls + p.lc > 0.0))
        throw std::invalid_argument("Tube: section lengths must be non-negative with positive total");
}

// ─── Getters ─────────────────────────────────────────────────────────────────

TubeParams Tube::parameters() const noexcept
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
    return m_u_ast[id];
}

double Tube::getK(Stiffness s) const noexcept
{
    const double csf = crossSectionFactor();
    switch (s)
    {
    case Stiffness::Bending:
        return m_E * pi_64 * csf;
    case Stiffness::Torsion:
        return m_G * pi_32 * csf;
    }
    return 0.0; // unreachable — satisfies compiler
}

} // namespace ctr
