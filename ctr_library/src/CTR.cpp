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

    // θᵢ(0) = αᵢ − βᵢ·u_iz(0) (twist wind-up over the straight transmission),
    // re-referenced so θ₁ ≡ 0 with α₁(0) absorbed into the base quaternion.
    m_theta_0 = {0.0, m_q[4UL] - b[1UL] * uz_0[1UL] - alpha1_0, m_q[5UL] - b[2UL] * uz_0[2UL] - alpha1_0};

    math::euler2Quaternion(0.0, alpha1_0, 0.0, m_h_0);
}

// ─── ODE integration (one forward shot) ──────────────────────────────────────

bvp_type CTR::residual(const bvp_type &initGuess, ShotMode mode)
{
    using namespace StateIdx;

    reset(initGuess);

    const bool record = (mode == ShotMode::Full);
    if (record)
    {
        m_y.clear();
        m_s.clear();
    }

    const auto &S = m_segment.getTransitionPoints();
    const blaze::StaticVector<double, NUM_TUBES> &distEnd = m_segment.getDistalEnds();

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

    const double ds = m_ds;
    const std::size_t n_seg = S.size() - 1UL;

    // uz₂/uz₃ at the respective tubes' distal ends. Distal ends are transition
    // points by construction, so they coincide with segment boundaries and are
    // captured inside the loop below for free. Initialize from s = 0 to cover
    // the degenerate fully-retracted case (distal end clamped to 0).
    m_uzDistal = {y_0[UZ_2], y_0[UZ_3]};

    constexpr double kBoundaryTol = 1.0e-7; // matches Segment's merge tolerance

    for (std::size_t seg_idx = 0UL; seg_idx < n_seg; ++seg_idx)
    {
        const double s_start = S[seg_idx];
        const double s_end = S[seg_idx + 1UL];

        m_stateEquations.setEquationParameters(m_segment, seg_idx, m_wf);

        // Fixed-step march: nFull steps of ds plus one explicit remainder step.
        const double len = s_end - s_start;
        const auto nFull = static_cast<std::size_t>(std::floor(len / ds + 1.0e-9));

        double s = s_start;
        if (record)
        {
            m_y.push_back(y_0);
            m_s.push_back(s);
        }
        for (std::size_t step = 0UL; step < nFull; ++step)
        {
            stepper.do_step(std::ref(m_stateEquations), y_0, s, ds);
            s += ds;
            if (record)
            {
                m_y.push_back(y_0);
                m_s.push_back(s);
            }
        }
        if (const double rem = s_end - s; rem > 1.0e-12)
        {
            stepper.do_step(std::ref(m_stateEquations), y_0, s, rem);
            if (record)
            {
                m_y.push_back(y_0);
                m_s.push_back(s_end);
            }
        }

        // Capture uz at tube 2/3 distal ends as we pass them.
        if (std::fabs(s_end - distEnd[1UL]) < kBoundaryTol)
            m_uzDistal[0UL] = y_0[UZ_2];
        if (std::fabs(s_end - distEnd[2UL]) < kBoundaryTol)
            m_uzDistal[1UL] = y_0[UZ_3];
    }

    m_finalState = y_0;

    // Distal boundary conditions, expressed as a NON-DIMENSIONAL residue: every
    // component is a curvature [1/m] (moment rows scaled by 1/EI₁, the torsion
    // row by 1/GJ₁), so the single L∞ tolerance weighs all rows equally and the
    // BVP Jacobian is well-scaled.
    //
    // Unloaded fast path: with no distal load the transverse proximal moment
    // vanishes exactly, so the moment rows are replaced by the equivalent
    // trivial constraints x₀ = x₁ = 0 (same solution set, exactly linear —
    // Newton satisfies them in one step and their Jacobian columns cost no
    // integrations).
    if (unloaded())
        return {initGuess[0UL], initGuess[1UL], y_0[UZ_1], m_uzDistal[0UL], m_uzDistal[1UL]};

    blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1;
    math::getSO3(blaze::subvector<QUAT_W, 4UL>(y_0), R1);

    const blaze::StaticVector<double, 3UL> distalMoment = blaze::trans(R1) * m_wm;

    return {(y_0[MB_X] - distalMoment[0UL]) / m_EI1, (y_0[MB_Y] - distalMoment[1UL]) / m_EI1,
            y_0[UZ_1] - (blaze::trans(ODESystem::kE3) * distalMoment) / m_GJ1, m_uzDistal[0UL], m_uzDistal[1UL]};
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

    // Unloaded fast path: the moment rows of the residual are exactly x₀ and
    // x₁, so their Jacobian columns are unit vectors — two of the five
    // integrations vanish.
    std::size_t start = 0UL;
    if (unloaded())
    {
        blaze::column(jac_bvp, 0UL) = blaze::StaticVector<double, BVP_DIM>{1.0, 0.0, 0.0, 0.0, 0.0};
        blaze::column(jac_bvp, 1UL) = blaze::StaticVector<double, BVP_DIM>{0.0, 1.0, 0.0, 0.0, 0.0};
        start = 2UL;
    }

    for (std::size_t iter = start; iter < BVP_DIM; ++iter)
    {
        initGuessPerturbed[iter] += scaled[iter];
        residuePerturbed = residual(initGuessPerturbed);
        blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / scaled[iter];
        initGuessPerturbed[iter] = initGuess[iter];
    }

    return jac_bvp;
}

Mat<3UL, 6UL> CTR::kinematicJacobian(const bvp_type &xStar)
{
    // Finite-difference steps. Small steps are viable ONLY because integration
    // is deterministic fixed-step: truncation error cancels between nominal
    // and perturbed shots (β-perturbations move segment boundaries, adding
    // O(ds⁴) ≈ 1e-12 discretization jitter — far below the step sizes below).
    constexpr double dX = 1.0e-6;     // shooting variables [1/m]
    constexpr double dBeta = 1.0e-5;  // [m]
    constexpr double dAlpha = 1.0e-4; // [rad]

    // Base shot at (q, x*): residual F₀ and tip r₀.
    const bvp_type F0 = residual(xStar);
    const blaze::StaticVector<double, 3UL> tip0 = tipPosition();

    // ∂F/∂x (5×5) and ∂r/∂x (3×5) — one lean shot per shooting variable.
    // Unloaded fast path: the moment rows of F are exactly x₀/x₁, so their
    // Fx columns are unit vectors and (because the corresponding rows of
    // Fx⁻¹Fq are then zero) the matching Rx columns never contribute — both
    // perturbation shots are skipped.
    Mat<BVP_DIM, BVP_DIM> Fx;
    Mat<3UL, BVP_DIM> Rx;
    std::size_t startX = 0UL;
    if (unloaded())
    {
        blaze::column(Fx, 0UL) = blaze::StaticVector<double, BVP_DIM>{1.0, 0.0, 0.0, 0.0, 0.0};
        blaze::column(Fx, 1UL) = blaze::StaticVector<double, BVP_DIM>{0.0, 1.0, 0.0, 0.0, 0.0};
        blaze::column(Rx, 0UL) = 0.0;
        blaze::column(Rx, 1UL) = 0.0;
        startX = 2UL;
    }

    bvp_type xp(xStar);
    for (std::size_t j = startX; j < BVP_DIM; ++j)
    {
        xp[j] += dX;
        const bvp_type Fp = residual(xp);
        blaze::column(Fx, j) = (Fp - F0) / dX;
        blaze::column(Rx, j) = (tipPosition() - tip0) / dX;
        xp[j] = xStar[j];
    }

    // ∂F/∂q (5×6) and ∂r/∂q (3×6) — one lean shot per joint, at FIXED x*.
    Mat<BVP_DIM, 6UL> Fq;
    Mat<3UL, 6UL> Rq;

    const blaze::StaticVector<double, 6UL> q_original(m_q);
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

    for (std::size_t i = 0UL; i < 6UL; ++i)
    {
        const double step = (i < NUM_TUBES) ? dBeta : dAlpha;
        m_q[i] = q_original[i] + step;

        // Angular DOFs do not alter segment transition points — but the FIRST
        // angular iteration (i == NUM_TUBES) must still recalculate once to
        // clear the β₃ perturbation left by iteration i = NUM_TUBES-1;
        // otherwise every α column is contaminated by (dBeta/dAlpha) times the
        // β₃ segment response.
        if (i <= NUM_TUBES)
            m_segment.recalculateSegments(m_tubes, betaView());

        const bvp_type Fp = residual(xStar);
        blaze::column(Fq, i) = (Fp - F0) / step;
        blaze::column(Rq, i) = (tipPosition() - tip0) / step;
        m_q[i] = q_original[i];
    }

    // Implicit function theorem on F(x, q) = const:
    //   dx/dq = −(∂F/∂x)⁻¹ ∂F/∂q  ⇒  dtip/dq = ∂r/∂q − (∂r/∂x)(∂F/∂x)⁻¹(∂F/∂q).
    const Mat<BVP_DIM, 6UL> Y = detail::invertJacobian(Fx) * Fq;
    return Rq - Rx * Y;
    // ScopeExit destructor runs here — m_q and m_segment are restored.
}

// ─── Actuation ────────────────────────────────────────────────────────────────

FKResult CTR::actuate(const blaze::StaticVector<double, 6UL> &q, bvp_type &initGuess)
{
    setConfiguration(q);
    m_segment.recalculateSegments(m_tubes, betaView());
    ensureSolver();
    ShootingProblem problem(*this);
    const FKResult result = m_solver->solve(initGuess, problem);

    // One Full shot at the returned shooting vector: the recorded backbone
    // trajectory (shapes, states) then always corresponds exactly to the
    // result the caller sees, regardless of which perturbed shot a solver
    // evaluated last. Costs one integration (~5-10% of a warm solve).
    std::ignore = residual(initGuess, ShotMode::Full);
    return result;
}

// ─── Inverse kinematics ───────────────────────────────────────────────────────

IKResult CTR::solveIK(const blaze::StaticVector<double, 3UL> &target, double posTol, bvp_type &initGuess,
                      const IKOptions &opts)
{
    detail::sanitizeBVPGuess(initGuess);

    // ── Establish a valid FK state at the current configuration ──────────────
    FKResult fk = actuate(m_q, initGuess);
    if (!fk)
    {
        initGuess = 0.0; // one clean cold restart
        fk = actuate(m_q, initGuess);
        if (!fk)
            return {.converged = false,
                    .positionError = blaze::norm(target - tipPosition()),
                    .iterations = 0UL,
                    .lastBVPStatus = fk.status,
                    .q = m_q};
    }

    // Accepted state: configuration, shooting vector, tip, error.
    blaze::StaticVector<double, 6UL> q_acc = m_q;
    bvp_type x_acc = initGuess;
    blaze::StaticVector<double, 3UL> tip_acc = tipPosition();
    blaze::StaticVector<double, 3UL> e = target - tip_acc;
    double err_acc = blaze::norm(e);

    // Tube geometry for the telescoping joint limits.
    constexpr double Clr = 5.0E-3; // clearance between tube ends [m]
    const blaze::StaticVector<double, NUM_TUBES> L = {m_tubes[0UL].getTubeLength(), m_tubes[1UL].getTubeLength(),
                                                      m_tubes[2UL].getTubeLength()};
    const blaze::StaticVector<double, NUM_TUBES> ls = {m_tubes[0UL].getStraightLen(), m_tubes[1UL].getStraightLen(),
                                                       m_tubes[2UL].getStraightLen()};

    // Levenberg-Marquardt damping state (dimensionally scaled inside the loop).
    double lambda = opts.dampingSeed;
    double nu = 2.0;

    std::size_t iter = 0UL;
    std::size_t stallCount = 0UL;

    while ((err_acc > posTol) && (iter < opts.maxIterations))
    {
        ++iter;

        // Re-establish the accepted configuration (a rejected trial leaves m_q
        // at the trial point) — cheap: no BVP solve, x_acc is already converged.
        setConfiguration(q_acc);
        m_segment.recalculateSegments(m_tubes, betaView());

        // ── Total kinematic Jacobian on the equilibrium manifold ─────────────
        Mat<3UL, 6UL> J = kinematicJacobian(x_acc);
        if (!blaze::isfinite(J))
        {
            // Rare: re-solve at the accepted configuration and retry once.
            detail::sanitizeBVPGuess(x_acc);
            fk = actuate(q_acc, x_acc);
            J = kinematicJacobian(x_acc);
            if (!fk || !blaze::isfinite(J))
                break;
        }

        // ── Damped least-squares step in the CAP-NORMALIZED joint space ──────
        // β [m] and α [rad] are incommensurable; treating them Euclidean lets
        // the large-magnitude β columns dominate the step direction, and a
        // uniform cap then starves the α components. Scaling each joint by its
        // own per-iteration cap (q̃ᵢ = qᵢ/sᵢ) makes "one unit" mean "one full
        // step" for every joint: J̃ = J·S, solve (J̃J̃ᵀ + μI)h = e, Δq̃ = J̃ᵀh,
        // and a uniform ‖Δq̃‖∞ ≤ 1 cap is fair across joint types.
        e = target - tip_acc;

        blaze::StaticVector<double, 6UL> S;
        for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
            S[i] = opts.maxBetaStep;
        for (std::size_t i = NUM_TUBES; i < 6UL; ++i)
            S[i] = opts.maxAlphaStep;

        Mat<3UL, 6UL> Js(J);
        for (std::size_t i = 0UL; i < 6UL; ++i)
            blaze::column(Js, i) *= S[i];

        const Mat<3UL, 3UL> JJt = Js * blaze::trans(Js);
        const double diagScale = blaze::max(blaze::abs(blaze::diagonal(JJt)));
        const double mu = lambda * diagScale;

        blaze::StaticVector<double, 3UL> h;
        Mat<3UL, 3UL> Adamped(JJt);
        for (std::size_t i = 0UL; i < 3UL; ++i)
            Adamped(i, i) += mu;
        blaze::solve(blaze::declsym(Adamped), h, e);

        blaze::StaticVector<double, 6UL> dqScaled = blaze::trans(Js) * h;

        // Uniform cap in scaled space: no component exceeds its own step cap.
        const double dqInf = blaze::linfNorm(dqScaled);
        if (dqInf > 1.0)
            dqScaled *= 1.0 / dqInf;

        blaze::StaticVector<double, 6UL> dq = S * dqScaled; // element-wise back-scaling

        // ── Trial configuration: telescoping β limits + α wrap ────────────────
        blaze::StaticVector<double, 6UL> q_trial = q_acc + dq;
        {
            const auto &b = q_acc; // limits evaluated at the accepted β
            blaze::StaticVector<double, NUM_TUBES> betaMin, betaMax;
            betaMin[0UL] = std::max({-ls[0UL], L[1UL] + b[1UL] - L[0UL], L[2UL] + b[2UL] - L[0UL]});
            betaMin[1UL] = std::max({-ls[1UL], b[0UL] + Clr, L[2UL] + b[2UL] - L[1UL]});
            betaMin[2UL] = std::max(-ls[2UL], b[1UL] + Clr);

            betaMax[0UL] = b[1UL] - Clr;
            betaMax[1UL] = std::min(b[2UL] - Clr, L[0UL] + b[0UL] - L[1UL]);
            betaMax[2UL] = std::min(L[1UL] + b[1UL] - L[2UL], L[0UL] + b[0UL] - L[2UL]);

            for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
                q_trial[i] = std::clamp(q_trial[i], betaMin[i], betaMax[i]);
        }
        // The wrap is an exact 2π congruence — FK-invariant, so the linear
        // model below keeps using the unwrapped dq.
        const blaze::StaticVector<double, 6UL> dq_actual = q_trial - q_acc;
        blaze::subvector<3UL, NUM_TUBES>(q_trial) =
            blaze::map(blaze::subvector<3UL, NUM_TUBES>(q_trial), [](double a) { return math::wrapToPi(a); });

        // ── Trial evaluation via warm-started BVP solve ───────────────────────
        bvp_type x_trial = x_acc;
        const FKResult fkTrial = actuate(q_trial, x_trial);

        bool accepted = false;
        if (fkTrial)
        {
            const blaze::StaticVector<double, 3UL> tip_trial = tipPosition();
            const double err_trial = blaze::norm(target - tip_trial);

            if (err_trial < err_acc)
            {
                // Gain ratio for the damping update (Madsen-Nielsen-Tingleff).
                const double predicted = blaze::sqrNorm(e) - blaze::sqrNorm(e - J * dq_actual);
                const double rho =
                    (blaze::sqrNorm(e) - err_trial * err_trial) / std::max(predicted, 1.0e-300);

                accepted = true;
                q_acc = q_trial;
                x_acc = x_trial;
                tip_acc = tip_trial;
                fk = fkTrial;

                if (err_acc - err_trial < 1.0e-10)
                    ++stallCount;
                else
                    stallCount = 0UL;
                err_acc = err_trial;

                lambda *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * rho - 1.0, 3.0));
                lambda = std::max(lambda, 1.0e-12);
                nu = 2.0;
            }
        }

        if (!accepted)
        {
            // Rejected (worse error or failed BVP): increase damping, shrink steps.
            lambda *= nu;
            nu *= 2.0;
            if (lambda > 1.0e8)
                break; // model exhausted — no productive step exists
        }

        if (stallCount >= 3UL)
            break; // three consecutive negligible improvements
    }

    // ── Exit invariant: object state corresponds to the accepted result ──────
    initGuess = x_acc;
    fk = actuate(q_acc, initGuess);

    const double finalError = blaze::norm(target - tipPosition());
    return {.converged = (finalError <= posTol),
            .positionError = finalError,
            .iterations = iter,
            .lastBVPStatus = fk.status,
            .q = m_q};
}

// ─── Shape access ─────────────────────────────────────────────────────────────

blaze::StaticVector<double, 3UL> CTR::tipPosition() const
{
    return blaze::subvector<StateIdx::POS_X, 3UL>(m_finalState);
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
