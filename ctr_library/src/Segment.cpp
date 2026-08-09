#include "ctr/Segment.hpp"
#include <algorithm>
#include <cmath>

namespace ctr
{

// ─── Constructors ─────────────────────────────────────────────────────────────

Segment::Segment(const std::array<Tube, NUM_TUBES> &tubes, const blaze::StaticVector<double, NUM_TUBES> &beta)
{
    recalculateSegments(tubes, beta);
}

// ─── Private: recalculate (called by CTR via friend) ──────────────────────────

void Segment::recalculateSegments(const std::array<Tube, NUM_TUBES> &tubes,
                                  const blaze::StaticVector<double, NUM_TUBES> &beta)
{
    m_S.clear();
    m_S.reserve(2UL * NUM_TUBES + 1UL);
    m_S.emplace_back(0.0);

    for (std::size_t i = 0; i < NUM_TUBES; ++i)
    {
        // Clamp to s >= 0: a tube's curved section (or, degenerately, the whole
        // tube) may be retracted into the actuation unit; the integrated
        // backbone always starts at the robot base s = 0.
        m_dist_end[i] = std::max(0.0, tubes[i].getTubeLength() + beta[i]);
        m_len_curv[i] = std::max(0.0, m_dist_end[i] - tubes[i].getCurvLen());

        m_S.emplace_back(m_len_curv[i]);
        m_S.emplace_back(m_dist_end[i]);
    }

    static constexpr double kTol = 1.0E-7;
    auto nearlyEqual = [](double a, double b) { return std::fabs(a - b) < kTol; };

    std::sort(m_S.begin(), m_S.end());
    m_S.erase(std::unique(m_S.begin(), m_S.end(), nearlyEqual), m_S.end());

    const std::size_t nSeg = m_S.size() - 1UL;

    m_EI.resize(NUM_TUBES, nSeg, false);
    m_GJ.resize(NUM_TUBES, nSeg, false);
    m_U_x.resize(NUM_TUBES, nSeg, false);
    m_U_y.resize(NUM_TUBES, nSeg, false);

    // Zero all matrices before selectively filling.
    m_EI = 0.0;
    m_GJ = 0.0;
    m_U_x = 0.0;
    m_U_y = 0.0;

    auto sBegin = m_S.begin();

    for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
    {
        auto it_b = std::lower_bound(sBegin, m_S.end(), m_len_curv[i] - kTol);
        auto it_c = std::lower_bound(sBegin, m_S.end(), m_dist_end[i] - kTol);

        const std::size_t b = static_cast<std::size_t>(std::distance(sBegin, it_b));
        const std::size_t c = static_cast<std::size_t>(std::distance(sBegin, it_c));

        // Fill bending and torsional stiffness for all segments where tube i is present.
        blaze::submatrix(m_EI, i, 0UL, 1UL, c) = tubes[i].getK(Stiffness::Bending);
        blaze::submatrix(m_GJ, i, 0UL, 1UL, c) = tubes[i].getK(Stiffness::Torsion);

        // Fill pre-curvatures only over the curved portion.
        if (b != c)
        {
            const std::size_t span = c - b;
            blaze::submatrix(m_U_x, i, b, 1UL, span) = tubes[i].get_u_ast(0UL);
            blaze::submatrix(m_U_y, i, b, 1UL, span) = tubes[i].get_u_ast(1UL);
        }
    }
}

// ─── Getters ─────────────────────────────────────────────────────────────────

const std::vector<double> &Segment::getTransitionPoints() const noexcept
{
    return m_S;
}

const blaze::StaticVector<double, NUM_TUBES> &Segment::getDistalEnds() const noexcept
{
    return m_dist_end;
}

const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &Segment::get_EI() const noexcept
{
    return m_EI;
}

const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &Segment::get_GJ() const noexcept
{
    return m_GJ;
}

const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &Segment::get_U_x() const noexcept
{
    return m_U_x;
}

const blaze::HybridMatrix<double, NUM_TUBES, MAX_SEGMENTS, blaze::columnMajor> &Segment::get_U_y() const noexcept
{
    return m_U_y;
}

} // namespace ctr
