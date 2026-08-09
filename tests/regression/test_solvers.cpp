// Cross-validation of the BVP solver suite and externally-loaded cases.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestRobot.hpp"

#include <array>
#include <cmath>
#include <vector>

using namespace ctr;
using Catch::Matchers::WithinAbs;

namespace
{

std::vector<blaze::StaticVector<double, 6UL>> configGrid()
{
    std::vector<blaze::StaticVector<double, 6UL>> grid;
    const std::array<double, 3UL> alphas2 = {0.0, math::deg2Rad(60.0), math::deg2Rad(170.0)};
    const std::array<double, 3UL> alphas3 = {0.0, math::deg2Rad(-45.0), math::deg2Rad(120.0)};

    for (const double a2 : alphas2)
        for (const double a3 : alphas3)
        {
            blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
            q[4UL] = a2;
            q[5UL] = a3;
            grid.push_back(q);
        }

    // One configuration with different retractions.
    blaze::StaticVector<double, 6UL> q = {
        -100.0e-3, -85.0e-3, -70.0e-3, 0.0, math::deg2Rad(90.0), math::deg2Rad(-120.0)};
    grid.push_back(q);
    return grid;
}

constexpr std::array<RootFindingMethod, 4UL> kAllMethods = {
    RootFindingMethod::ModifiedNewtonRaphson, RootFindingMethod::LevenbergMarquardt, RootFindingMethod::PowellDogLeg,
    RootFindingMethod::Broyden};

} // namespace

TEST_CASE("All solvers agree on the converged shooting solution", "[regression][solvers]")
{
    const auto grid = configGrid();

    for (const auto &q : grid)
    {
        CAPTURE(q[3UL], q[4UL], q[5UL], q[0UL]);

        // Reference: default solver.
        CTR reference = testing::makeReferenceRobot();
        bvp_type guessRef{};
        const FKResult fkRef = reference.actuate(q, guessRef);
        REQUIRE(fkRef);
        const auto tipRef = reference.tipPosition();

        for (const RootFindingMethod method : kAllMethods)
        {
            CAPTURE(static_cast<int>(method));
            CTR robot = testing::makeReferenceRobot(method);
            bvp_type guess{};
            const FKResult fk = robot.actuate(q, guess);
            REQUIRE(fk);
            CHECK(fk.residual <= testing::kBVPTol);

            // Physical agreement: tips must match across solvers.
            CHECK(testing::maxAbsDiff(robot.tipPosition(), tipRef) < 5.0e-6);
            // Shooting-variable agreement (different solvers stop at different
            // points within the residual tolerance).
            CHECK(blaze::linfNorm(guess - guessRef) < 1.0e-3);
        }
    }
}

TEST_CASE("Externally loaded robot converges (distal force)", "[regression][solvers][loaded]")
{
    for (const RootFindingMethod method :
         {RootFindingMethod::ModifiedNewtonRaphson, RootFindingMethod::LevenbergMarquardt})
    {
        CAPTURE(static_cast<int>(method));

        CTR robot = testing::makeReferenceRobot(method);
        bvp_type guess{};
        REQUIRE(robot.actuate(testing::kHomeConfig, guess));
        const auto tipUnloaded = robot.tipPosition();

        robot.setDistalForce({0.0, 0.02, 0.0}); // 20 mN lateral tip load
        bvp_type guessLoaded{};
        const FKResult fk = robot.actuate(testing::kHomeConfig, guessLoaded);
        REQUIRE(fk);
        CHECK(fk.residual <= testing::kBVPTol);

        // The load must visibly deflect the tip...
        CHECK(testing::maxAbsDiff(robot.tipPosition(), tipUnloaded) > 1.0e-4);
        // ...and the tip must remain physically reachable.
        CHECK(blaze::norm(robot.tipPosition()) <= 0.13 + 1e-9);
    }
}

TEST_CASE("Loaded robot with distal moment converges", "[regression][solvers][loaded]")
{
    CTR robot = testing::makeReferenceRobot();
    robot.setDistalMoment({0.0, 2.0e-3, 0.0}); // 2 mN·m distal moment

    bvp_type guess{};
    const FKResult fk = robot.actuate(testing::kHomeConfig, guess);
    REQUIRE(fk);
    CHECK(fk.residual <= testing::kBVPTol);
    CHECK(blaze::isfinite(robot.tipPosition()));
}

TEST_CASE("Hard cold-start configuration: honest failure and Broyden escape", "[regression][solvers]")
{
    // alpha = [-30, 100, -160] deg is a documented hard cold start: the
    // Newton-type iterations descend into a non-root local minimum of ||f||^2
    // (residual ~0.097). Broyden's quasi-Newton path escapes it, which is why
    // the default solver's fallback cascade tries a zero-restarted Broyden.
    const blaze::StaticVector<double, 6UL> qHard = {
        -0.12, -0.1, -0.08, math::deg2Rad(-30.0), math::deg2Rad(100.0), math::deg2Rad(-160.0)};

    SECTION("Pure Newton-type solvers fail honestly (no fake convergence, finite residual)")
    {
        for (const RootFindingMethod method : {RootFindingMethod::LevenbergMarquardt, RootFindingMethod::PowellDogLeg})
        {
            CAPTURE(static_cast<int>(method));
            CTR robot = testing::makeReferenceRobot(method);
            bvp_type guess{};
            const FKResult fk = robot.actuate(qHard, guess);
            CHECK_FALSE(fk);
            CHECK(std::isfinite(fk.residual));
            CHECK(fk.residual > testing::kBVPTol);
        }
    }

    SECTION("The default solver converges via its Broyden fallback stage")
    {
        CTR robot = testing::makeReferenceRobot(RootFindingMethod::ModifiedNewtonRaphson);
        bvp_type guess{};
        const FKResult fk = robot.actuate(qHard, guess);
        REQUIRE(fk);
        CHECK(fk.residual <= testing::kBVPTol);
    }

    SECTION("Broyden converges cold")
    {
        CTR robot = testing::makeReferenceRobot(RootFindingMethod::Broyden);
        bvp_type guess{};
        const FKResult fk = robot.actuate(qHard, guess);
        REQUIRE(fk);
        CHECK(fk.residual <= testing::kBVPTol);
    }

    SECTION("Warm-start continuation reaches the configuration")
    {
        CTR robot = testing::makeReferenceRobot();
        bvp_type guess{};
        bool ok = true;
        for (int step = 1; step <= 4 && ok; ++step)
        {
            blaze::StaticVector<double, 6UL> qi = qHard;
            const double t = step / 4.0;
            qi[3UL] *= t;
            qi[4UL] *= t;
            qi[5UL] *= t;
            ok = static_cast<bool>(robot.actuate(qi, guess));
        }
        CHECK(ok);
    }
}

TEST_CASE("Unloaded fast path matches the full 5-unknown solve", "[regression][solvers]")
{
    // The unloaded fast path replaces the moment rows by the trivial
    // constraints x0 = x1 = 0 (3 shooting unknowns effectively). Forcing the
    // full 5-unknown path with a numerically negligible distal force must
    // give the same physical solution.
    for (const auto &q : configGrid())
    {
        CAPTURE(q[3UL], q[4UL], q[5UL]);

        CTR fast = testing::makeReferenceRobot();
        bvp_type gFast{};
        REQUIRE(fast.actuate(q, gFast));

        CTR full = testing::makeReferenceRobot();
        full.setDistalForce({0.0, 0.0, 1.0e-30}); // disables the fast path, physically nil
        bvp_type gFull{};
        REQUIRE(full.actuate(q, gFull));

        // Both formulations stop within the BVP tolerance of the same root;
        // the residual stopping-point difference shows up at the ~1e-8 level.
        CHECK(testing::maxAbsDiff(fast.tipPosition(), full.tipPosition()) < 5.0e-8);
        CHECK(blaze::linfNorm(gFast - gFull) < 1.0e-5);
    }
}

TEST_CASE("Re-actuating at the converged guess is a fixed point", "[regression][solvers]")
{
    // Guards the line-search bookkeeping: the trajectory recorded by the solver
    // must correspond exactly to the shooting vector it returns.
    CTR robot = testing::makeReferenceRobot();
    blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
    q[4UL] = math::deg2Rad(60.0);
    q[5UL] = math::deg2Rad(120.0);

    bvp_type guess{};
    REQUIRE(robot.actuate(q, guess));
    const auto tip1 = robot.tipPosition();
    const bvp_type guess1 = guess;

    const FKResult fk2 = robot.actuate(q, guess);
    REQUIRE(fk2);
    CHECK(fk2.iterations == 0UL); // already converged — no Newton step needed
    CHECK(testing::maxAbsDiff(robot.tipPosition(), tip1) < 1.0e-12);
    CHECK(blaze::linfNorm(guess - guess1) < 1.0e-12);
}
