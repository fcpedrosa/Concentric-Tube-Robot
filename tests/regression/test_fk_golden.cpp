// Golden-value FK regression tests.
//
// The golden data (golden_fk_data.hpp) was captured after the correctness-fix
// and RK4 fixed-step phases, from cold zero-guess solves with the default
// solver. These tests pin the FK numerics: any change that moves a tip beyond
// the tolerances below must be an intentional, documented numerical change
// with a justified golden re-capture.

#include <catch2/catch_test_macros.hpp>

#include "TestRobot.hpp"
#include "golden_fk_data.hpp"

#include <algorithm>

using namespace ctr;

namespace
{
// Pure refactors must reproduce goldens essentially exactly; intentional
// numerical changes (e.g. BVP rescaling) may relax these with justification.
constexpr double kTipTol = 1.0e-11;  // [m]
constexpr double kQuatTol = 1.0e-9;  // quaternion components (sign-normalized)
constexpr double kXStarTol = 1.0e-6; // shooting variables (mixed units until rescaling)
} // namespace

TEST_CASE("Golden FK: configuration grid incl. loaded cases", "[regression][golden]")
{
    for (std::size_t i = 0; i < testing::golden::kNumCases; ++i)
    {
        const auto &g = testing::golden::kCases[i];
        CAPTURE(i, g.q[0UL], g.q[3UL], g.q[4UL], g.q[5UL], g.force[1UL], g.moment[1UL]);

        CTR robot = testing::makeReferenceRobot();
        robot.setDistalForce(g.force);
        robot.setDistalMoment(g.moment);

        bvp_type guess{};
        const FKResult fk = robot.actuate(g.q, guess);
        REQUIRE(fk);
        CHECK(fk.residual <= testing::kBVPTol);

        CHECK(blaze::linfNorm(robot.tipPosition() - g.tip) < kTipTol);

        // Quaternion comparison up to global sign (q and -q are the same rotation).
        const auto states = robot.states();
        const auto &yEnd = states[states.size() - 1UL];
        blaze::StaticVector<double, 4UL> quat = {yEnd[StateIdx::QUAT_W], yEnd[StateIdx::QUAT_X],
                                                 yEnd[StateIdx::QUAT_Y], yEnd[StateIdx::QUAT_Z]};
        const double dPlus = blaze::linfNorm(quat - g.quat);
        const double dMinus = blaze::linfNorm(quat + g.quat);
        CHECK(std::min(dPlus, dMinus) < kQuatTol);

        CHECK(blaze::linfNorm(guess - g.xStar) < kXStarTol);
    }
}
