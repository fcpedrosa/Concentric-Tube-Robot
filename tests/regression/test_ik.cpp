// Inverse-kinematics regression tests: Jacobian correctness, FK -> IK -> FK
// round trips, convergence budgets, and honest failure semantics.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestRobot.hpp"

#include <vector>

using namespace ctr;

namespace
{
constexpr double kPosTol = 5.0e-4; // 0.5 mm, matching the demo
}

TEST_CASE("Total kinematic Jacobian matches central differences with full BVP re-solve", "[regression][ik]")
{
    // Reference: dtip/dq by central differences where each perturbed
    // configuration is fully RE-SOLVED (warm-started), i.e. the true total
    // derivative along the equilibrium manifold — exactly what the IFT-based
    // kinematicJacobian claims to compute.
    std::vector<blaze::StaticVector<double, 6UL>> configs;
    {
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        configs.push_back(q);
        q[4UL] = math::deg2Rad(60.0);
        q[5UL] = math::deg2Rad(120.0);
        configs.push_back(q);
        q = {-100.0e-3, -85.0e-3, -70.0e-3, math::deg2Rad(30.0), math::deg2Rad(-60.0), math::deg2Rad(90.0)};
        configs.push_back(q);
    }

    for (const auto &q0 : configs)
    {
        CAPTURE(q0[0UL], q0[3UL], q0[4UL], q0[5UL]);

        CTR robot = testing::makeReferenceRobot();
        bvp_type guess{};
        REQUIRE(robot.actuate(q0, guess));

        const Mat<3UL, 6UL> J = robot.kinematicJacobian(guess);
        REQUIRE(blaze::isfinite(J));

        for (std::size_t i = 0UL; i < 6UL; ++i)
        {
            const double delta = (i < NUM_TUBES) ? 1.0e-4 : 1.0e-3;

            blaze::StaticVector<double, 6UL> qp = q0, qm = q0;
            qp[i] += delta;
            qm[i] -= delta;

            bvp_type gp = guess, gm = guess;
            CTR rp = robot, rm = robot;
            REQUIRE(rp.actuate(qp, gp));
            REQUIRE(rm.actuate(qm, gm));

            const blaze::StaticVector<double, 3UL> col = (rp.tipPosition() - rm.tipPosition()) / (2.0 * delta);
            const double colNorm = blaze::norm(col);
            const double err = blaze::norm(blaze::column(J, i) - col);

            CAPTURE(i, colNorm, err);
            // Relative agreement where the column is significant; absolute
            // agreement for near-null columns. The reference itself carries
            // solver stopping-point noise of ~5e-4 absolute (tips resolved to
            // the BVP tolerance, divided by the 2δ baseline), so the gate
            // cannot be tighter than ~1e-3 relative on O(0.5) columns.
            if (colNorm > 1.0e-3)
                CHECK(err / colNorm < 2.5e-3);
            else
                CHECK(err < 1.0e-5);
        }
    }
}

TEST_CASE("IK: FK -> IK -> FK round trips", "[regression][ik]")
{
    // Sample reachable targets by running FK at known configurations, then
    // demand solveIK reaches each from the home configuration.
    std::vector<blaze::StaticVector<double, 6UL>> targetConfigs;
    {
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        q[4UL] = math::deg2Rad(30.0);
        targetConfigs.push_back(q);
        q[4UL] = math::deg2Rad(75.0);
        q[5UL] = math::deg2Rad(-45.0);
        targetConfigs.push_back(q);
        q = testing::kHomeConfig;
        q[0UL] += 8.0e-3; // retract tube 1 slightly less
        q[3UL] = math::deg2Rad(20.0);
        targetConfigs.push_back(q);
        q = {-110.0e-3, -95.0e-3, -78.0e-3, 0.0, math::deg2Rad(-60.0), math::deg2Rad(60.0)};
        targetConfigs.push_back(q);
        q = {-120.0e-3, -100.0e-3, -80.0e-3, math::deg2Rad(45.0), math::deg2Rad(45.0), math::deg2Rad(100.0)};
        targetConfigs.push_back(q);
    }

    std::size_t totalIterations = 0UL;
    for (const auto &qStar : targetConfigs)
    {
        CAPTURE(qStar[0UL], qStar[3UL], qStar[4UL], qStar[5UL]);

        // Target tip from FK.
        CTR probe = testing::makeReferenceRobot();
        bvp_type guessProbe{};
        REQUIRE(probe.actuate(qStar, guessProbe));
        const auto target = probe.tipPosition();

        // Solve IK from home.
        CTR robot = testing::makeReferenceRobot();
        bvp_type guess{};
        REQUIRE(robot.actuate(testing::kHomeConfig, guess));

        const IKResult ik = robot.solveIK(target, kPosTol, guess);
        CAPTURE(ik.positionError, ik.iterations);
        CHECK(ik.converged);
        CHECK(ik.positionError <= kPosTol);
        CHECK(ik.iterations <= 30UL);

        // The returned state is self-consistent.
        CHECK(testing::maxAbsDiff(robot.tipPosition(), target) <= kPosTol * 2.0);
        CHECK(robot.configuration() == ik.q);

        totalIterations += ik.iterations;
    }
    CAPTURE(totalIterations);
    CHECK(totalIterations <= 100UL);
}

TEST_CASE("IK: reachable target requiring a large telescoping excursion", "[regression][ik]")
{
    // The round trips above all keep β within 10 mm of home, so they never
    // exercise a long β travel — nor the telescoping limits that bound it.
    // This target (the demo's) needs 25/30/45 mm of β and drives the iterate
    // onto the tube-clearance face, where a naive per-coordinate clamp against
    // the neighbours' pre-step β deadlocks: β₁ cannot advance until β₂ has and
    // vice versa. That crawl left the tip 1.2 mm out after 500 iterations.
    const blaze::StaticVector<double, 6UL> qStar = {
        -95.0e-3, -70.0e-3, -35.0e-3, 0.0, math::deg2Rad(60.0), math::deg2Rad(-75.0)};

    CTR probe = testing::makeReferenceRobot();
    bvp_type guessProbe{};
    REQUIRE(probe.actuate(qStar, guessProbe)); // the target is reachable by construction
    const auto target = probe.tipPosition();

    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));

    const IKResult ik = robot.solveIK(target, kPosTol, guess);
    CAPTURE(ik.positionError, ik.iterations);
    CHECK(ik.converged);
    CHECK(ik.positionError <= kPosTol);
    CHECK(testing::maxAbsDiff(robot.tipPosition(), target) <= kPosTol * 2.0);

    // A real budget, not just "under the cap": the deadlock burned every
    // iteration it was given, so a generous-but-finite bound still catches a
    // convergence-rate regression.
    CHECK(ik.iterations <= 60UL);

    // The returned configuration must satisfy the telescoping limits. The old
    // per-coordinate clamp could violate them outright — when two tubes moved
    // in the same step, each was bounded by the other's stale value, so the
    // pair could close past the clearance.
    constexpr double kClr = 5.0e-3; // must match the clearance used by solveIK
    constexpr double kSlop = 1.0e-9;
    const auto &q = ik.q;
    CAPTURE(q[0UL], q[1UL], q[2UL]);
    CHECK(q[1UL] - q[0UL] >= kClr - kSlop);
    CHECK(q[2UL] - q[1UL] >= kClr - kSlop);
    for (std::size_t i = 0UL; i < NUM_TUBES; ++i)
        CHECK(q[i] >= -robot.tubes()[i].getStraightLen() - kSlop);
    // Inner tubes must reach at least as far as the ones outside them.
    const auto &tubes = robot.tubes();
    CHECK(tubes[0UL].getTubeLength() + q[0UL] >= tubes[1UL].getTubeLength() + q[1UL] - kSlop);
    CHECK(tubes[1UL].getTubeLength() + q[1UL] >= tubes[2UL].getTubeLength() + q[2UL] - kSlop);
}

TEST_CASE("IK: unreachable target reports failure honestly", "[regression][ik]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));

    // Far outside the workspace (backbone is at most 0.13 m long).
    const blaze::StaticVector<double, 3UL> unreachable = {0.0, 0.0, 0.5};

    const IKResult ik = robot.solveIK(unreachable, kPosTol, guess);
    CHECK_FALSE(ik.converged);
    CHECK(std::isfinite(ik.positionError));
    CHECK(ik.positionError > 0.3); // cannot get closer than ~0.37 m
    CHECK(blaze::isfinite(robot.tipPosition()));

    // The object is left in a consistent, converged FK state.
    CHECK(ik.lastBVPStatus == SolverStatus::Converged);
}

TEST_CASE("IK: repeated solves warm-start and remain stable", "[regression][ik]")
{
    // Track a short sequence of nearby targets (2 mm apart) — the tracking
    // use-case that motivates warm starts.
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));

    blaze::StaticVector<double, 3UL> target = robot.tipPosition();
    for (std::size_t i = 1UL; i <= 5UL; ++i)
    {
        target[1UL] += 2.0e-3;
        const IKResult ik = robot.solveIK(target, kPosTol, guess);
        CAPTURE(i, ik.positionError, ik.iterations);
        CHECK(ik.converged);
        CHECK(ik.iterations <= 10UL);
    }
}
