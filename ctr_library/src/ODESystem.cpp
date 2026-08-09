#include "ctr/ODESystem.hpp"

#include <cmath>

namespace ctr
{

// default constructor
ODESystem::ODESystem() : m_u_ast_x(0.0), m_u_ast_y(0.0), m_GJ(0.0), m_EIoverGJ(0.0), m_EIux(0.0), m_EIuy(0.0), m_f(0.0)
{
}

void ODESystem::setEquationParameters(const Segment &seg, std::size_t segIdx,
                                      const blaze::StaticVector<double, 3UL> &force)
{
    const auto &EI = seg.get_EI();
    const auto &GJ = seg.get_GJ();
    const auto &U_x = seg.get_U_x();
    const auto &U_y = seg.get_U_y();

    double sumEI = 0.0;
    for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
    {
        const double ei = EI(i, segIdx);
        const double gj = GJ(i, segIdx);
        m_u_ast_x[i] = U_x(i, segIdx);
        m_u_ast_y[i] = U_y(i, segIdx);
        m_GJ[i] = gj;
        m_EIoverGJ[i] = (gj != 0.0) ? ei / gj : 0.0;
        m_EIux[i] = ei * m_u_ast_x[i];
        m_EIuy[i] = ei * m_u_ast_y[i];
        sumEI += ei;
    }
    // At least tube 1 is present in every segment, so sumEI > 0.
    m_invSumEI = 1.0 / sumEI;

    m_f = force;
    m_hasForce = (blaze::sqrNorm(force) != 0.0);
}

// Cosserat-model ODE right-hand side for a three-tube CTR (Rucker, Jones &
// Webster, IEEE T-RO 2010), written in scalar form: evaluated four times per
// integration step, so no matrix temporaries are constructed here.
void ODESystem::operator()(const state_type &y, state_type &dyds, const double /*s*/) const noexcept
{
    using namespace StateIdx;

    // Relative twist angles of tubes 2 and 3 w.r.t. tube 1.
    const double c2 = std::cos(y[THETA_2]), s2 = std::sin(y[THETA_2]);
    const double c3 = std::cos(y[THETA_3]), s3 = std::sin(y[THETA_3]);

    // Bending moment of the assembly, body frame: mb = [y0, y1, Σ GJᵢ·uzᵢ].
    const double mbx = y[MB_X];
    const double mby = y[MB_Y];
    const double mbz = m_GJ[0UL] * y[UZ_1] + m_GJ[1UL] * y[UZ_2] + m_GJ[2UL] * y[UZ_3];

    // Curvature of tube 1 (xy): u1_xy = (mb + Σ Rz(θᵢ)·Kᵢ·u*ᵢ)_xy / ΣEI,
    // with the Rz products expanded in scalars (Kᵢ·u*ᵢ has zero z-component).
    const double u1x = m_invSumEI * (mbx + m_EIux[0UL] + (c2 * m_EIux[1UL] - s2 * m_EIuy[1UL]) +
                                     (c3 * m_EIux[2UL] - s3 * m_EIuy[2UL]));
    const double u1y = m_invSumEI * (mby + m_EIuy[0UL] + (s2 * m_EIux[1UL] + c2 * m_EIuy[1UL]) +
                                     (s3 * m_EIux[2UL] + c3 * m_EIuy[2UL]));
    const double u1z = y[UZ_1];

    // Curvatures of tubes 2 and 3 (xy): uᵢ = Rz(θᵢ)ᵀ·u1 + θ'ᵢ·e3.
    const double u2x = c2 * u1x + s2 * u1y;
    const double u2y = -s2 * u1x + c2 * u1y;
    const double u3x = c3 * u1x + s3 * u1y;
    const double u3y = -s3 * u1x + c3 * u1y;

    // Torsional curvature rates u'zᵢ = (EIᵢ/GJᵢ)(uₓᵢ·u*ᵧᵢ − uᵧᵢ·u*ₓᵢ) and twist
    // rates θ'ᵢ = uzᵢ − uz₁ (identically 0 for tube 1 and for absent tubes,
    // whose m_EIoverGJ and m_GJ are 0).
    dyds[UZ_1] = m_EIoverGJ[0UL] * (u1x * m_u_ast_y[0UL] - u1y * m_u_ast_x[0UL]);
    dyds[UZ_2] = m_EIoverGJ[1UL] * (u2x * m_u_ast_y[1UL] - u2y * m_u_ast_x[1UL]);
    dyds[UZ_3] = m_EIoverGJ[2UL] * (u3x * m_u_ast_y[2UL] - u3y * m_u_ast_x[2UL]);
    dyds[THETA_1] = 0.0;
    dyds[THETA_2] = (m_GJ[1UL] != 0.0) ? (y[UZ_2] - y[UZ_1]) : 0.0;
    dyds[THETA_3] = (m_GJ[2UL] != 0.0) ? (y[UZ_3] - y[UZ_1]) : 0.0;

    // Moment rate (xy): mb' = −u1 × mb − ê₃·R₁ᵀ·f  (force term only when loaded).
    dyds[MB_X] = -(u1y * mbz - u1z * mby);
    dyds[MB_Y] = -(u1z * mbx - u1x * mbz);

    // Quaternion h = [w, x, y, z] of the tube-1 body frame.
    const double hw = y[QUAT_W], hx = y[QUAT_X], hy = y[QUAT_Y], hz = y[QUAT_Z];
    // Same self-normalizing scale as math::getSO3 (absorbs quaternion drift).
    const double scale = 2.0 / (hw * hw + hx * hx + hy * hy + hz * hz);

    if (m_hasForce)
    {
        // (R₁ᵀ f) needs the full rotation; ê₃·v = e₃ × v = [−v_y, v_x, 0].
        const double r00 = 1.0 + scale * (-hy * hy - hz * hz);
        const double r10 = scale * (hx * hy + hz * hw);
        const double r20 = scale * (hx * hz - hy * hw);
        const double r01 = scale * (hx * hy - hz * hw);
        const double r11 = 1.0 + scale * (-hx * hx - hz * hz);
        const double r21 = scale * (hy * hz + hx * hw);
        const double r02 = scale * (hx * hz + hy * hw);
        const double r12 = scale * (hy * hz - hx * hw);
        const double r22 = 1.0 + scale * (-hx * hx - hy * hy);

        const double vx = r00 * m_f[0UL] + r10 * m_f[1UL] + r20 * m_f[2UL]; // (R₁ᵀ f)_x
        const double vy = r01 * m_f[0UL] + r11 * m_f[1UL] + r21 * m_f[2UL]; // (R₁ᵀ f)_y
        dyds[MB_X] -= -vy;
        dyds[MB_Y] -= vx;

        // r' = R₁·e₃ (third column).
        dyds[POS_X] = r02;
        dyds[POS_Y] = r12;
        dyds[POS_Z] = r22;
    }
    else
    {
        // Only the third column of R₁ is needed.
        dyds[POS_X] = scale * (hx * hz + hy * hw);
        dyds[POS_Y] = scale * (hy * hz - hx * hw);
        dyds[POS_Z] = 1.0 + scale * (-hx * hx - hy * hy);
    }

    // Quaternion rate: h' = 0.5 · h ⊗ [0, u1].
    dyds[QUAT_W] = 0.5 * (-u1x * hx - u1y * hy - u1z * hz);
    dyds[QUAT_X] = 0.5 * (u1x * hw + u1z * hy - u1y * hz);
    dyds[QUAT_Y] = 0.5 * (u1y * hw - u1z * hx + u1x * hz);
    dyds[QUAT_Z] = 0.5 * (u1z * hw + u1y * hx - u1x * hy);
}

} // namespace ctr
