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
            // agreement for near-null columns.
            if (colNorm > 1.0e-3)
                CHECK(err / colNorm < 1.0e-3);
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
