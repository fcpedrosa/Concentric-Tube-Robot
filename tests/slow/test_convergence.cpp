// Integration convergence study: verifies the fixed-step RK4 integrator's
// observed order and that the default step ds = 1 mm delivers tip accuracy
// far below the IK position tolerance.

#include <catch2/catch_test_macros.hpp>

#include "TestRobot.hpp"

#include <array>
#include <cmath>
#include <vector>

using namespace ctr;

namespace
{

// Warm-started across step sizes: the study measures integration accuracy of
// one physical solution, not cold-start solver robustness (covered elsewhere).
blaze::StaticVector<double, 3UL> tipAt(const blaze::StaticVector<double, 6UL> &q, double ds, bvp_type &guess)
{
    CTR robot = testing::makeReferenceRobot();
    robot.setIntegrationStep(ds);
    REQUIRE(robot.actuate(q, guess));
    return robot.tipPosition();
}

std::vector<blaze::StaticVector<double, 6UL>> studyGrid()
{
    std::vector<blaze::StaticVector<double, 6UL>> grid;
    for (const double a2 : {0.0, math::deg2Rad(60.0), math::deg2Rad(150.0)})
        for (const double a3 : {math::deg2Rad(-45.0), math::deg2Rad(120.0)})
        {
            blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
            q[4UL] = a2;
            q[5UL] = a3;
            grid.push_back(q);
        }
    return grid;
}

} // namespace

TEST_CASE("Fixed-step RK4: observed order and default-step accuracy", "[slow][convergence]")
{
    const auto grid = studyGrid();

    for (const auto &q : grid)
    {
        CAPTURE(q[4UL], q[5UL]);

        bvp_type guess{};
        const auto reference = tipAt(q, 0.25e-3, guess);

        const std::array<double, 3UL> steps = {4.0e-3, 2.0e-3, 1.0e-3};
        std::array<double, 3UL> errors{};
        for (std::size_t i = 0; i < steps.size(); ++i)
            errors[i] = blaze::linfNorm(tipAt(q, steps[i], guess) - reference);

        CAPTURE(errors[0], errors[1], errors[2]);

        // Default-step accuracy: 1 mm step must land within 1e-7 m of the
        // refined reference (posTol is 5e-4 m — 3+ orders of margin).
        // Measured: ~3e-11 m, i.e. RK4 truncation error is negligible here.
        CHECK(errors[2] < 1.0e-7);

        // Refining the step must not make things worse (allowing noise-floor
        // slack: near ~1e-10 m the BVP tolerance dominates, so a clean
        // asymptotic order-4 slope is not observable — absolute accuracy is
        // the meaningful gate).
        CHECK(errors[0] < 1.0e-6);
        CHECK(errors[1] <= errors[0] * 2.0);
        CHECK(errors[2] <= errors[1] * 2.0);
    }
}
