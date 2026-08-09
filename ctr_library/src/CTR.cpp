#include "CTR.hpp"
#include "BVPSolver.hpp"

#include <algorithm>
#include <cmath>

namespace ctr
{

// ─── Constructor & copy/move ─────────────────────────────────────────────────

CTR::CTR(const std::array<std::shared_ptr<Tube>, NUM_TUBES> &Tb, const blaze::StaticVector<double, 6UL> &q, double Tol,
         mathOp::rootFindingMethod method)
    : m_Tubes(Tb), m_q(q), m_theta_0(0.0), m_h_0{1.0, 0.0, 0.0, 0.0}, m_wf(0.0), m_wm(0.0), m_accuracy(Tol),
      m_method(method), m_segment(rawTubes(), blaze::subvector<0UL, NUM_TUBES>(q)), m_stateEquations(),
      m_solver(makeBVPSolver(method))
{
    m_theta_0 = {0.0, q[4UL] - q[3UL], q[5UL] - q[4UL]};
    m_y.reserve(1000UL);
    m_s.reserve(1000UL);
}

CTR::CTR(const CTR &rhs)
    : m_Tubes(rhs.m_Tubes), m_q(rhs.m_q), m_theta_0(rhs.m_theta_0), m_h_0(rhs.m_h_0), m_wf(rhs.m_wf), m_wm(rhs.m_wm),
      m_accuracy(rhs.m_accuracy), m_method(rhs.m_method),
      m_segment(rhs.rawTubes(), blaze::subvector<0UL, NUM_TUBES>(rhs.m_q)), m_stateEquations(), m_y(rhs.m_y),
      m_s(rhs.m_s), m_solver(makeBVPSolver(rhs.m_method))
{
}

CTR &CTR::operator=(const CTR &rhs)
{
    if (this != &rhs)
    {
        m_Tubes = rhs.m_Tubes;
        m_q = rhs.m_q;
        m_theta_0 = rhs.m_theta_0;
        m_h_0 = rhs.m_h_0;
        m_wf = rhs.m_wf;
        m_wm = rhs.m_wm;
        m_accuracy = rhs.m_accuracy;
        m_method = rhs.m_method;
        m_segment = Segment(rhs.rawTubes(), blaze::subvector<0UL, NUM_TUBES>(rhs.m_q));
        m_stateEquations = ODESystem{};
        m_y = rhs.m_y;
        m_s = rhs.m_s;
        m_solver = makeBVPSolver(rhs.m_method);
    }
    return *this;
}

// ─── ODE reset ───────────────────────────────────────────────────────────────

void CTR::reset(const bvp_type &initGuess)
{
    using namespace StateIdx;
    blaze::StaticVector<double, NUM_TUBES> uz_0 = {initGuess[UZ_1], initGuess[UZ_2], initGuess[UZ_3]};
    const auto b = beta();
    const double alpha1_0 = m_q[3UL] - b[0UL] * uz_0[0UL];

    m_y.clear();
    m_y.reserve(1000UL);
    m_s.clear();
    m_s.reserve(1000UL);

    m_theta_0 = {0.0, m_q[4UL] - b[1UL] * uz_0[1UL] - alpha1_0, m_q[5UL] - b[2UL] * uz_0[2UL] - alpha1_0};

    mathOp::euler2Quaternion(0.0, alpha1_0, 0.0, m_h_0);
}

// ─── ODE integration ─────────────────────────────────────────────────────────

bvp_type CTR::ODESolver(const bvp_type &initGuess)
{
    using namespace StateIdx;

    reset(initGuess);

    const SegmentData seg = m_segment.getParameters();
    const auto &[EI, GJ, U_x, U_y, S] = seg;

    boost::numeric::odeint::adaptive_adams_bashforth_moulton<8UL, State, double, State, double, BlazeBVPAlgebra,
                                                             boost::numeric::odeint::default_operations,
                                                             boost::numeric::odeint::initially_resizer>
        abm8_stepper;

    // Initial conditions
    state_type y_0;
    y_0[MB_X] = initGuess[0UL];
    y_0[MB_Y] = initGuess[1UL];
    y_0[UZ_1] = initGuess[2UL];
    y_0[UZ_2] = initGuess[3UL];
    y_0[UZ_3] = initGuess[4UL];
    y_0[THETA_1] = m_theta_0[0UL];
    y_0[THETA_2] = m_theta_0[1UL];
    y_0[THETA_3] = m_theta_0[2UL];
    y_0[POS_X] = 0.0;
    y_0[POS_Y] = 0.0;
    y_0[POS_Z] = 0.0;
    y_0[QUAT_W] = m_h_0[0UL];
    y_0[QUAT_X] = m_h_0[1UL];
    y_0[QUAT_Y] = m_h_0[2UL];
    y_0[QUAT_Z] = m_h_0[3UL];

    auto observer = [this](const state_type &y, double s)
    {
        m_y.push_back(y);
        m_s.push_back(s);
    };

    constexpr double ds = 1.0E-3;
    const std::size_t n_seg = S.size() - 1UL;

    for (std::size_t seg_idx = 0UL; seg_idx < n_seg; ++seg_idx)
    {
        const double s_start = S[seg_idx];
        const double s_end = S[seg_idx + 1UL];

        m_stateEquations.setEquationParameters(blaze::column(U_x, seg_idx), blaze::column(U_y, seg_idx),
                                               blaze::column(EI, seg_idx), blaze::column(GJ, seg_idx), m_wf);

        boost::numeric::odeint::integrate_adaptive(abm8_stepper, m_stateEquations, y_0, s_start, s_end, ds, observer);
    }

    // Distal boundary conditions
    blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1;
    mathOp::getSO3(blaze::subvector<QUAT_W, 4UL>(y_0), R1);

    const blaze::StaticVector<double, 3UL> distalMoment = blaze::trans(R1) * m_wm;

    bvp_type Residue = {y_0[MB_X] - distalMoment[0UL], y_0[MB_Y] - distalMoment[1UL],
                        GJ(0UL, 0UL) * y_0[UZ_1] - blaze::trans(ODESystem::kE3) * distalMoment, 0.0, 0.0};

    const blaze::StaticVector<double, NUM_TUBES> distEnd = m_segment.getDistalEnds();

    auto computeResidue = [&](double distalEnd, std::size_t index) -> void
    {
        auto itt = std::lower_bound(m_s.begin(), m_s.end(), distalEnd - 1.0E-7);
        const auto id = static_cast<std::size_t>(std::distance(m_s.begin(), itt));
        Residue[2UL + index] = m_y[id][UZ_1 + index];
    };

    computeResidue(distEnd[1UL], 1UL);
    computeResidue(distEnd[2UL], 2UL);

    return Residue;
}

// ─── Jacobians ────────────────────────────────────────────────────────────────

blaze::StaticMatrix<double, BVP_DIM, BVP_DIM> CTR::jac_BVP(const bvp_type &initGuess, const bvp_type &residue)
{
    blaze::StaticMatrix<double, BVP_DIM, BVP_DIM, blaze::columnMajor> jac_bvp;

    bvp_type initGuessPerturbed(initGuess), residuePerturbed, scaled(initGuess);
    constexpr double incr_scale = 1.0E-7;
    constexpr double incr_floor = 1.0E-9;

    scaled *= incr_scale;
    scaled = blaze::generate(
        BVP_DIM, [&](std::size_t idx)
        { return (std::fabs(scaled[idx]) > incr_floor) ? scaled[idx] : std::copysign(incr_floor, initGuess[idx]); });

    for (std::size_t iter = 0UL; iter < BVP_DIM; ++iter)
    {
        initGuessPerturbed[iter] += scaled[iter];
        residuePerturbed = ODESolver(initGuessPerturbed);
        blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / scaled[iter];
        initGuessPerturbed[iter] = initGuess[iter];
    }

    return jac_bvp;
}

blaze::StaticMatrix<double, 3UL, 6UL> CTR::jacobian(const bvp_type &initGuess,
                                                    const blaze::StaticVector<double, 3UL> &tipPos)
{
    blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> jac;

    const blaze::StaticVector<double, 6UL> q_original(m_q);

    constexpr double incr_scale = 1.0E-3;
    constexpr double incr_floor = 5.0E-4;

    blaze::StaticVector<double, 6UL> q_scaled(m_q);
    q_scaled *= incr_scale;
    q_scaled = blaze::generate(6UL,
                               [&](std::size_t idx)
                               {
                                   return (std::fabs(q_scaled[idx]) > incr_floor)
                                              ? q_scaled[idx]
                                              : std::copysign(incr_floor, q_original[idx]);
                               });

    // RAII guard: restores m_q and m_segment on any exit path (normal or exception).
    auto doRestore = [&]() noexcept
    {
        m_q = q_original;
        m_segment.recalculateSegments(rawTubes(), beta());
    };
    struct ScopeExit
    {
        decltype(doRestore) &fn;
        ~ScopeExit() noexcept { fn(); }
    } guard{doRestore};

    blaze::StaticVector<double, 6UL> q_perturbed(m_q);
    for (std::size_t iter = 0UL; iter <= 5UL; ++iter)
    {
        q_perturbed[iter] += q_scaled[iter];
        m_q = q_perturbed;

        // Angular DOFs (iter >= NUM_TUBES) do not alter segment transition points.
        if (iter < NUM_TUBES)
            m_segment.recalculateSegments(rawTubes(), beta());

        std::ignore = ODESolver(initGuess);
        blaze::column(jac, iter) = (getTipPos() - tipPos) / q_scaled[iter];
        q_perturbed[iter] = q_original[iter];
    }

    return jac;
    // ScopeExit destructor runs here — m_q and m_segment are restored.
}

// ─── Actuation ────────────────────────────────────────────────────────────────

bool CTR::actuate_CTR(bvp_type &initGuess, const blaze::StaticVector<double, 6UL> &q_input)
{
    setConfiguration(q_input);
    m_segment.recalculateSegments(rawTubes(), beta());
    return m_solver->solve(initGuess, *this);
}

// ─── Position control ─────────────────────────────────────────────────────────

bool CTR::posCTRL(bvp_type &initGuess, const blaze::StaticVector<double, 3UL> &target, double posTol)
{
    double minError = 1.0E3;
    bool status;
    blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> J;
    blaze::StaticMatrix<double, 6UL, 3UL, blaze::columnMajor> J_inv;

    detail::sanitizeBVPGuess(initGuess);

    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>> Kp, Kd, Ki;
    blaze::diagonal(Kp) = 1.000;
    blaze::diagonal(Ki) = 0.050;
    blaze::diagonal(Kd) = 0.001;

    blaze::StaticVector<double, 6UL> dqdt, q_min(m_q), q(m_q);
    bvp_type initGuessMin(initGuess);
    status = actuate_CTR(initGuess, q);

    if (!status)
        return status;

    blaze::StaticVector<double, 3UL> x_CTR = getTipPos();
    blaze::StaticVector<double, 3UL> tipError = target - x_CTR;
    blaze::StaticVector<double, 3UL> last_tipError = tipError;
    blaze::StaticVector<double, 3UL> d_tipError{0.0, 0.0, 0.0};
    blaze::StaticVector<double, 3UL> int_tipError{0.0, 0.0, 0.0};
    static constexpr double integralClamp = 5.0;

    double dist2Tgt = blaze::norm(tipError);

    if (dist2Tgt < minError)
    {
        minError = dist2Tgt;
        q_min = q;
        if (dist2Tgt <= posTol)
            return status;
    }

    blaze::StaticVector<double, 6UL> f{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    constexpr double Clr = 5.0E-3;

    const blaze::StaticVector<double, NUM_TUBES> L = {m_Tubes[0UL]->getTubeLength(), m_Tubes[1UL]->getTubeLength(),
                                                      m_Tubes[2UL]->getTubeLength()};

    const blaze::StaticVector<double, NUM_TUBES> ls = {m_Tubes[0UL]->getStraightLen(), m_Tubes[1UL]->getStraightLen(),
                                                       m_Tubes[2UL]->getStraightLen()};

    blaze::StaticVector<double, NUM_TUBES> betaMax, betaMin;

    std::size_t N_itr = 0UL;
    constexpr std::size_t maxIter = 500UL;
    constexpr double ke = 4.0;

    while ((dist2Tgt > posTol) && (N_itr < maxIter))
    {
        N_itr++;

        J = jacobian(initGuess, x_CTR);

        std::size_t isfinite_ctr = 0UL;
        while (!blaze::isfinite(J))
        {
            ++isfinite_ctr;
            if (isfinite_ctr > 200UL)
                break;
            initGuess *= 0.750;
            detail::sanitizeBVPGuess(initGuess);
            std::ignore = ODESolver(initGuess);
            x_CTR = getTipPos();
            J = jacobian(initGuess, x_CTR);
        }

        J_inv = mathOp::pInv(J);

        const auto b = beta();
        betaMin[0UL] = std::max({-ls[0UL], L[1UL] + b[1UL] - L[0UL], L[2UL] + b[2UL] - L[0UL]});
        betaMin[1UL] = std::max({-ls[1UL], b[0UL] + Clr, L[2UL] + b[2UL] - L[1UL]});
        betaMin[2UL] = std::max(-ls[2UL], b[1UL] + Clr);

        betaMax[0UL] = b[1UL] - Clr;
        betaMax[1UL] = std::min(b[2UL] - Clr, L[0UL] + b[0UL] - L[1UL]);
        betaMax[2UL] = std::min(L[1UL] + b[1UL] - L[2UL], L[0UL] + b[0UL] - L[2UL]);

        constexpr double spanEps = 1.0E-6;
        blaze::StaticVector<double, NUM_TUBES> span = betaMax - betaMin;
        for (std::size_t i = 0; i < NUM_TUBES; ++i)
        {
            const double magnitude = std::max(std::abs(span[i]), spanEps);
            span[i] = std::copysign(magnitude, span[i]);
        }
        const auto normalizedOffset = (betaMax + betaMin - 2.0 * b) / span;
        blaze::subvector<0UL, NUM_TUBES>(f) =
            blaze::pow(blaze::abs(normalizedOffset), ke) * blaze::sign(b - (betaMax + betaMin) * 0.5);

        const auto taskSpaceCommand = Kp * tipError + Kd * d_tipError + Ki * int_tipError;
        dqdt = J_inv * taskSpaceCommand;

        auto rescale_dqdt = [&]() noexcept -> void
        {
            const auto b2 = beta();
            for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
            {
                if (b2[i] + dqdt[i] > betaMax[i])
                    dqdt[i] = (betaMax[i] - b2[i]) * 0.5;
                if (b2[i] + dqdt[i] < betaMin[i])
                    dqdt[i] = (betaMin[i] - b2[i]) * 0.5;
            }
        };

        rescale_dqdt();

        q += dqdt;
        blaze::subvector<3UL, NUM_TUBES>(q) =
            blaze::map(blaze::subvector<3UL, NUM_TUBES>(q), [](double theta) { return mathOp::congruentAngle(theta); });

        status = actuate_CTR(initGuess, q);

        if (!status)
        {
            initGuess *= 0.75;
            detail::sanitizeBVPGuess(initGuess);
            status = actuate_CTR(initGuess, q);

            if (!status)
            {
                initGuess = initGuessMin;
                std::ignore = actuate_CTR(initGuess, q_min);
            }
        }

        x_CTR = getTipPos();
        tipError = target - x_CTR;

        int_tipError += tipError;
        int_tipError =
            blaze::map(int_tipError, [](double value) { return std::clamp(value, -integralClamp, integralClamp); });

        d_tipError = tipError - last_tipError;
        last_tipError = tipError;
        dist2Tgt = blaze::norm(tipError);

        if (dist2Tgt < minError)
        {
            minError = dist2Tgt;
            q_min = q;
            initGuessMin = initGuess;
        }

        if (blaze::linfNorm(dqdt) <= 1.0E-8)
        {
            initGuess = initGuessMin;
            return actuate_CTR(initGuess, q_min);
        }
    }

    initGuess = std::move(initGuessMin);
    return actuate_CTR(initGuess, q_min);
}

// ─── Getters ─────────────────────────────────────────────────────────────────

std::array<std::shared_ptr<Tube>, NUM_TUBES> CTR::getTubes() const
{
    return m_Tubes;
}

blaze::StaticVector<double, 3UL> CTR::getBeta() const
{
    return blaze::subvector<0UL, NUM_TUBES>(m_q);
}

blaze::StaticVector<double, 6UL> CTR::getConfiguration() const
{
    return m_q;
}

blaze::StaticVector<double, 3UL> CTR::getTipPos() const
{
    blaze::StaticVector<double, 3UL> pos;
    if (!m_y.empty())
        pos = blaze::subvector<StateIdx::POS_X, 3UL>(m_y.back());
    return pos;
}

blaze::StaticVector<double, NUM_TUBES> CTR::getDistalEnds() const
{
    return m_segment.getDistalEnds();
}

std::tuple<blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>,
           blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>,
           blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor>>
CTR::getTubeShapes() const
{
    blaze::HybridMatrix<double, 3UL, 1000UL, blaze::columnMajor> Tb_1(3UL, m_y.size());

    const blaze::StaticVector<double, NUM_TUBES> distalEnds = m_segment.getDistalEnds();

    for (std::size_t col = 0UL; col < Tb_1.columns(); ++col)
    {
        blaze::column(Tb_1, col) = blaze::subvector<StateIdx::POS_X, 3UL>(m_y[col]) * 1.0E3;
    }

    auto tubeEndIndex = [&](std::size_t tube_index) -> std::size_t
    {
        auto it = std::lower_bound(m_s.begin(), m_s.end(), distalEnds[tube_index] - 1.0E-7);
        return static_cast<std::size_t>(std::distance(m_s.begin(), it));
    };

    const std::size_t distalIndex_Tb2 = tubeEndIndex(1UL);
    const std::size_t distalIndex_Tb3 = tubeEndIndex(2UL);

    auto Tb_2 = blaze::submatrix(Tb_1, 0UL, 0UL, 3UL, distalIndex_Tb2 + 1UL);
    auto Tb_3 = blaze::submatrix(Tb_1, 0UL, 0UL, 3UL, distalIndex_Tb3 + 1UL);

    return {Tb_1, Tb_2, Tb_3};
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> CTR::getShape() const
{
    std::vector<double> r_x, r_y, r_z;
    r_x.reserve(m_y.size());
    r_y.reserve(m_y.size());
    r_z.reserve(m_y.size());

    for (const auto &el : m_y)
    {
        r_x.emplace_back(el[StateIdx::POS_X]);
        r_y.emplace_back(el[StateIdx::POS_Y]);
        r_z.emplace_back(el[StateIdx::POS_Z]);
    }

    return {std::move(r_x), std::move(r_y), std::move(r_z)};
}

std::span<const state_type> CTR::states() const noexcept
{
    return {m_y.data(), m_y.size()};
}

std::span<const double> CTR::arcLengthSamples() const noexcept
{
    return {m_s.data(), m_s.size()};
}

// ─── Setters ─────────────────────────────────────────────────────────────────

void CTR::setConfiguration(const blaze::StaticVector<double, 6UL> &q)
{
    m_q = q;
}

void CTR::setBVPMethod(mathOp::rootFindingMethod mthd)
{
    m_method = mthd;
    m_solver = makeBVPSolver(mthd);
}

void CTR::setDistalMoment(const blaze::StaticVector<double, 3UL> &moment)
{
    m_wm = moment;
}

void CTR::setDistalForce(const blaze::StaticVector<double, 3UL> &force)
{
    m_wf = force;
}

} // namespace ctr
