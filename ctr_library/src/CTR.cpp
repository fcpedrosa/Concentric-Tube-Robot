#include "ctr/CTR.hpp"
#include "ctr/detail/mathOperations.hpp"
#include "BVPSolver.hpp"

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <tuple>

namespace ctr
{

// ─── Constructor & special members ───────────────────────────────────────────

CTR::CTR(std::array<Tube, NUM_TUBES> tubes, const blaze::StaticVector<double, 6UL> &q, double bvpTolerance,
         RootFindingMethod method)
    : m_tubes(std::move(tubes)), m_q(q), m_theta_0(0.0), m_h_0{1.0, 0.0, 0.0, 0.0}, m_wf(0.0), m_wm(0.0),
      m_accuracy(bvpTolerance), m_method(method), m_segment(m_tubes, blaze::subvector<0UL, NUM_TUBES>(q)),
      m_stateEquations(), m_solver(makeBVPSolver(method))
{
    m_theta_0 = {0.0, q[4UL] - q[3UL], q[5UL] - q[4UL]};
    m_EI1 = m_tubes[0UL].getK(Stiffness::Bending);
    m_GJ1 = m_tubes[0UL].getK(Stiffness::Torsion);
    m_y.reserve(1000UL);
    m_s.reserve(1000UL);
}

CTR::CTR(const CTR &rhs)
    : m_tubes(rhs.m_tubes), m_q(rhs.m_q), m_theta_0(rhs.m_theta_0), m_h_0(rhs.m_h_0), m_wf(rhs.m_wf), m_wm(rhs.m_wm),
      m_accuracy(rhs.m_accuracy), m_method(rhs.m_method), m_segment(rhs.m_segment),
      m_stateEquations(rhs.m_stateEquations), m_y(rhs.m_y), m_s(rhs.m_s), m_solver(makeBVPSolver(rhs.m_method))
{
}

CTR &CTR::operator=(const CTR &rhs)
{
    if (this != &rhs)
    {
        m_tubes = rhs.m_tubes;
        m_q = rhs.m_q;
        m_theta_0 = rhs.m_theta_0;
        m_h_0 = rhs.m_h_0;
        m_wf = rhs.m_wf;
        m_wm = rhs.m_wm;
        m_accuracy = rhs.m_accuracy;
        m_method = rhs.m_method;
        m_segment = rhs.m_segment;
        m_stateEquations = rhs.m_stateEquations;
        m_y = rhs.m_y;
        m_s = rhs.m_s;
        m_solver = makeBVPSolver(rhs.m_method);
    }
    return *this;
}

// Moves are defaulted out-of-line (BVPSolver must be complete for ~unique_ptr).
// A moved-from CTR has a null m_solver; every entry point that needs the solver
// calls ensureSolver(), which lazily recreates it from the retained m_method.
CTR::CTR(CTR &&rhs) noexcept = default;
CTR &CTR::operator=(CTR &&rhs) noexcept = default;
CTR::~CTR() = default;

void CTR::ensureSolver()
{
    if (!m_solver) [[unlikely]]
        m_solver = makeBVPSolver(m_method);
}

// ─── ShootingProblem facade ───────────────────────────────────────────────────

bvp_type ShootingProblem::residual(const bvp_type &x)
{
    return m_ctr.residual(x);
}

Mat<BVP_DIM, BVP_DIM> ShootingProblem::jacobian(const bvp_type &x, const bvp_type &f0)
{
    return m_ctr.bvpJacobian(x, f0);
}

double ShootingProblem::tolerance() const noexcept
{
    return m_ctr.tolerance();
}

// ─── ODE reset ───────────────────────────────────────────────────────────────

void CTR::reset(const bvp_type &initGuess)
{
    using namespace StateIdx;
    const blaze::StaticVector<double, NUM_TUBES> uz_0 = {initGuess[UZ_1], initGuess[UZ_2], initGuess[UZ_3]};
    const auto b = betaView();
    const double alpha1_0 = m_q[3UL] - b[0UL] * uz_0[0UL];

    m_y.clear();
    m_s.clear();

    // θᵢ(0) = αᵢ − βᵢ·u_iz(0) (twist wind-up over the straight transmission),
    // re-referenced so θ₁ ≡ 0 with α₁(0) absorbed into the base quaternion.
    m_theta_0 = {0.0, m_q[4UL] - b[1UL] * uz_0[1UL] - alpha1_0, m_q[5UL] - b[2UL] * uz_0[2UL] - alpha1_0};

    math::euler2Quaternion(0.0, alpha1_0, 0.0, m_h_0);
}

// ─── ODE integration (one forward shot) ──────────────────────────────────────

bvp_type CTR::residual(const bvp_type &initGuess)
{
    using namespace StateIdx;

    reset(initGuess);

    const auto &EI = m_segment.get_EI();
    const auto &GJ = m_segment.get_GJ();
    const auto &U_x = m_segment.get_U_x();
    const auto &U_y = m_segment.get_U_y();
    const auto &S = m_segment.getTransitionPoints();

    // Classic RK4, driven by a hand-written fixed-step loop below.
    //
    // Deliberately NOT an error-controlled stepper: both finite-difference
    // Jacobians (bvpJacobian and the task-space Jacobian) rely on the
    // truncation error cancelling between nominal and perturbed shots, which
    // requires a deterministic step sequence. A one-step method also has no
    // multistep restart transient at the 5-6 material discontinuities per
    // backbone (the previous ABM8 spent most short segments inside its
    // initialization ramp — and its error estimate was silently discarded by
    // the integrate_adaptive overload resolution anyway).
    boost::numeric::odeint::runge_kutta4<state_type, double, state_type, double,
                                         boost::numeric::odeint::vector_space_algebra>
        stepper;

    // Initial conditions. The shooting vector is non-dimensionalized: its
    // moment components are curvatures mb/EI₁ [1/m]; the ODE state carries
    // physical moments [N·m].
    state_type y_0;
    y_0[MB_X] = initGuess[0UL] * m_EI1;
    y_0[MB_Y] = initGuess[1UL] * m_EI1;
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

    auto observe = [this](const state_type &y, double s)
    {
        m_y.push_back(y);
        m_s.push_back(s);
    };

    const double ds = m_ds;
    const std::size_t n_seg = S.size() - 1UL;

    for (std::size_t seg_idx = 0UL; seg_idx < n_seg; ++seg_idx)
    {
        const double s_start = S[seg_idx];
        const double s_end = S[seg_idx + 1UL];

        m_stateEquations.setEquationParameters(blaze::column(U_x, seg_idx), blaze::column(U_y, seg_idx),
                                               blaze::column(EI, seg_idx), blaze::column(GJ, seg_idx), m_wf);

        // Fixed-step march: nFull steps of ds plus one explicit remainder step.
        const double len = s_end - s_start;
        const auto nFull = static_cast<std::size_t>(std::floor(len / ds + 1.0e-9));

        double s = s_start;
        observe(y_0, s);
        for (std::size_t step = 0UL; step < nFull; ++step)
        {
            stepper.do_step(std::ref(m_stateEquations), y_0, s, ds);
            s += ds;
            observe(y_0, s);
        }
        if (const double rem = s_end - s; rem > 1.0e-12)
        {
            stepper.do_step(std::ref(m_stateEquations), y_0, s, rem);
            observe(y_0, s_end);
        }
    }

    // Distal boundary conditions, expressed as a NON-DIMENSIONAL residue: every
    // component is a curvature [1/m] (moment rows scaled by 1/EI₁, the torsion
    // row by 1/GJ₁), so the single L∞ tolerance weighs all rows equally and the
    // BVP Jacobian is well-scaled.
    blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1;
    math::getSO3(blaze::subvector<QUAT_W, 4UL>(y_0), R1);

    const blaze::StaticVector<double, 3UL> distalMoment = blaze::trans(R1) * m_wm;

    bvp_type Residue = {(y_0[MB_X] - distalMoment[0UL]) / m_EI1, (y_0[MB_Y] - distalMoment[1UL]) / m_EI1,
                        y_0[UZ_1] - (blaze::trans(ODESystem::kE3) * distalMoment) / m_GJ1, 0.0, 0.0};

    const blaze::StaticVector<double, NUM_TUBES> &distEnd = m_segment.getDistalEnds();

    auto computeResidue = [&](double distalEnd, std::size_t index) -> void
    {
        auto itt = std::lower_bound(m_s.begin(), m_s.end(), distalEnd - 1.0E-7);
        // Clamp: lower_bound may return end() if the distal end exceeds the last
        // sample by more than the tolerance (defensive — S.back() == max distal end).
        const auto id = std::min(static_cast<std::size_t>(std::distance(m_s.begin(), itt)), m_s.size() - 1UL);
        Residue[2UL + index] = m_y[id][UZ_1 + index];
    };

    computeResidue(distEnd[1UL], 1UL);
    computeResidue(distEnd[2UL], 2UL);

    return Residue;
}

// ─── Jacobians ────────────────────────────────────────────────────────────────

Mat<BVP_DIM, BVP_DIM> CTR::bvpJacobian(const bvp_type &initGuess, const bvp_type &residue)
{
    Mat<BVP_DIM, BVP_DIM> jac_bvp;

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
        residuePerturbed = residual(initGuessPerturbed);
        blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / scaled[iter];
        initGuessPerturbed[iter] = initGuess[iter];
    }

    return jac_bvp;
}

Mat<3UL, 6UL> CTR::kinematicJacobian(const bvp_type &initGuess, const blaze::StaticVector<double, 3UL> &tipPos)
{
    Mat<3UL, 6UL> jac;

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
        m_segment.recalculateSegments(m_tubes, betaView());
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
            m_segment.recalculateSegments(m_tubes, betaView());

        std::ignore = residual(initGuess);
        blaze::column(jac, iter) = (tipPosition() - tipPos) / q_scaled[iter];
        q_perturbed[iter] = q_original[iter];
    }

    return jac;
    // ScopeExit destructor runs here — m_q and m_segment are restored.
}

// ─── Actuation ────────────────────────────────────────────────────────────────

FKResult CTR::actuate(const blaze::StaticVector<double, 6UL> &q, bvp_type &initGuess)
{
    setConfiguration(q);
    m_segment.recalculateSegments(m_tubes, betaView());
    ensureSolver();
    ShootingProblem problem(*this);
    return m_solver->solve(initGuess, problem);
}

// ─── Inverse kinematics ───────────────────────────────────────────────────────

IKResult CTR::solveIK(const blaze::StaticVector<double, 3UL> &target, double posTol, bvp_type &initGuess)
{
    double minError = 1.0E3;
    Mat<3UL, 6UL> J;
    Mat<6UL, 3UL> J_inv;

    detail::sanitizeBVPGuess(initGuess);

    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor>> Kp, Kd, Ki;
    blaze::diagonal(Kp) = 1.000;
    blaze::diagonal(Ki) = 0.050;
    blaze::diagonal(Kd) = 0.001;

    blaze::StaticVector<double, 6UL> dqdt, q_min(m_q), q(m_q);
    bvp_type initGuessMin(initGuess);

    FKResult fk = actuate(q, initGuess);
    if (!fk)
        return {.converged = false,
                .positionError = blaze::norm(target - tipPosition()),
                .iterations = 0UL,
                .lastBVPStatus = fk.status,
                .q = m_q};

    blaze::StaticVector<double, 3UL> x_CTR = tipPosition();
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
            return {.converged = true,
                    .positionError = dist2Tgt,
                    .iterations = 0UL,
                    .lastBVPStatus = fk.status,
                    .q = m_q};
    }

    constexpr double Clr = 5.0E-3;

    const blaze::StaticVector<double, NUM_TUBES> L = {m_tubes[0UL].getTubeLength(), m_tubes[1UL].getTubeLength(),
                                                      m_tubes[2UL].getTubeLength()};

    const blaze::StaticVector<double, NUM_TUBES> ls = {m_tubes[0UL].getStraightLen(), m_tubes[1UL].getStraightLen(),
                                                       m_tubes[2UL].getStraightLen()};

    blaze::StaticVector<double, NUM_TUBES> betaMax, betaMin;

    std::size_t N_itr = 0UL;
    constexpr std::size_t maxIter = 500UL;

    while ((dist2Tgt > posTol) && (N_itr < maxIter))
    {
        N_itr++;

        J = kinematicJacobian(initGuess, x_CTR);

        std::size_t isfinite_ctr = 0UL;
        while (!blaze::isfinite(J))
        {
            ++isfinite_ctr;
            if (isfinite_ctr > 200UL)
                break;
            initGuess *= 0.750;
            detail::sanitizeBVPGuess(initGuess);
            std::ignore = residual(initGuess);
            x_CTR = tipPosition();
            J = kinematicJacobian(initGuess, x_CTR);
        }

        J_inv = math::pInv(J);

        const auto b = betaView();
        betaMin[0UL] = std::max({-ls[0UL], L[1UL] + b[1UL] - L[0UL], L[2UL] + b[2UL] - L[0UL]});
        betaMin[1UL] = std::max({-ls[1UL], b[0UL] + Clr, L[2UL] + b[2UL] - L[1UL]});
        betaMin[2UL] = std::max(-ls[2UL], b[1UL] + Clr);

        betaMax[0UL] = b[1UL] - Clr;
        betaMax[1UL] = std::min(b[2UL] - Clr, L[0UL] + b[0UL] - L[1UL]);
        betaMax[2UL] = std::min(L[1UL] + b[1UL] - L[2UL], L[0UL] + b[0UL] - L[2UL]);

        const auto taskSpaceCommand = Kp * tipError + Kd * d_tipError + Ki * int_tipError;
        dqdt = J_inv * taskSpaceCommand;

        auto rescale_dqdt = [&]() noexcept -> void
        {
            const auto b2 = betaView();
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
            blaze::map(blaze::subvector<3UL, NUM_TUBES>(q), [](double theta) { return math::wrapToPi(theta); });

        fk = actuate(q, initGuess);

        if (!fk)
        {
            initGuess *= 0.75;
            detail::sanitizeBVPGuess(initGuess);
            fk = actuate(q, initGuess);

            if (!fk)
            {
                initGuess = initGuessMin;
                std::ignore = actuate(q_min, initGuess);
            }
        }

        x_CTR = tipPosition();
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
            break; // stalled — exit to the common wrap-up below
    }

    // Re-actuate at the best configuration found so the object state (shape,
    // tip, m_q) corresponds exactly to the returned result.
    initGuess = std::move(initGuessMin);
    fk = actuate(q_min, initGuess);

    const double finalError = blaze::norm(target - tipPosition());
    return {.converged = (finalError <= posTol),
            .positionError = finalError,
            .iterations = N_itr,
            .lastBVPStatus = fk.status,
            .q = m_q};
}

// ─── Shape access ─────────────────────────────────────────────────────────────

blaze::StaticVector<double, 3UL> CTR::tipPosition() const
{
    blaze::StaticVector<double, 3UL> pos;
    if (!m_y.empty())
        pos = blaze::subvector<StateIdx::POS_X, 3UL>(m_y.back());
    return pos;
}

blaze::StaticVector<double, NUM_TUBES> CTR::distalEnds() const
{
    return m_segment.getDistalEnds();
}

CTR::Points CTR::backboneShape() const
{
    Points shape;
    shape.reserve(m_y.size());
    for (const auto &y : m_y)
        shape.emplace_back(blaze::subvector<StateIdx::POS_X, 3UL>(y));
    return shape;
}

std::array<CTR::Points, NUM_TUBES> CTR::tubeShapes() const
{
    Points tube1 = backboneShape();

    const blaze::StaticVector<double, NUM_TUBES> &distal = m_segment.getDistalEnds();

    auto tubeEndIndex = [&](std::size_t tube_index) -> std::size_t
    {
        auto it = std::lower_bound(m_s.begin(), m_s.end(), distal[tube_index] - 1.0E-7);
        const auto id = static_cast<std::size_t>(std::distance(m_s.begin(), it));
        return (m_s.empty()) ? 0UL : std::min(id, m_s.size() - 1UL);
    };

    Points tube2, tube3;
    if (!tube1.empty())
    {
        tube2.assign(tube1.begin(), tube1.begin() + static_cast<std::ptrdiff_t>(tubeEndIndex(1UL) + 1UL));
        tube3.assign(tube1.begin(), tube1.begin() + static_cast<std::ptrdiff_t>(tubeEndIndex(2UL) + 1UL));
    }

    return {std::move(tube1), std::move(tube2), std::move(tube3)};
}

// ─── Observers ────────────────────────────────────────────────────────────────

const std::array<Tube, NUM_TUBES> &CTR::tubes() const noexcept
{
    return m_tubes;
}

blaze::StaticVector<double, 6UL> CTR::configuration() const noexcept
{
    return m_q;
}

blaze::StaticVector<double, NUM_TUBES> CTR::beta() const noexcept
{
    return blaze::subvector<0UL, NUM_TUBES>(m_q);
}

std::span<const state_type> CTR::states() const noexcept
{
    return {m_y.data(), m_y.size()};
}

std::span<const double> CTR::arcLengthSamples() const noexcept
{
    return {m_s.data(), m_s.size()};
}

// ─── Modifiers ────────────────────────────────────────────────────────────────

void CTR::setConfiguration(const blaze::StaticVector<double, 6UL> &q)
{
    m_q = q;
}

void CTR::setBVPMethod(RootFindingMethod method)
{
    m_method = method;
    m_solver = makeBVPSolver(method);
}

void CTR::setDistalMoment(const blaze::StaticVector<double, 3UL> &moment)
{
    m_wm = moment;
}

void CTR::setDistalForce(const blaze::StaticVector<double, 3UL> &force)
{
    m_wf = force;
}

void CTR::setTube(std::size_t idx, Tube tube)
{
    m_tubes[idx] = std::move(tube);
    m_segment.recalculateSegments(m_tubes, betaView());
    m_EI1 = m_tubes[0UL].getK(Stiffness::Bending);
    m_GJ1 = m_tubes[0UL].getK(Stiffness::Torsion);
}

void CTR::setIntegrationStep(double ds)
{
    m_ds = std::clamp(ds, 1.0e-5, 1.0e-2);
}

} // namespace ctr
