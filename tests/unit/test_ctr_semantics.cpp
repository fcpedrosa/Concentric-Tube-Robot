#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestRobot.hpp"

#include <utility>

using namespace ctr;
using Catch::Matchers::WithinAbs;

TEST_CASE("FK at home converges and reports a sane tip", "[unit][ctr]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};

    const FKResult fk = robot.actuate(testing::kHomeConfig, guess);
    REQUIRE(fk);
    CHECK(fk.residual <= testing::kBVPTol);

    const auto tip = robot.tipPosition();
    CHECK(blaze::isfinite(tip));
    // Backbone length at home is 0.13 m; the tip cannot be farther than that.
    CHECK(blaze::norm(tip) <= 0.13 + 1e-9);
    CHECK(blaze::norm(tip) > 0.05);
}

TEST_CASE("Copies are independent", "[unit][ctr]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));
    const auto tipBefore = robot.tipPosition();

    CTR copy(robot);
    blaze::StaticVector<double, 6UL> q2 = testing::kHomeConfig;
    q2[3UL] = math::deg2Rad(45.0);
    bvp_type guess2{};
    REQUIRE(copy.actuate(q2, guess2));

    // The original must be untouched by actuating the copy.
    CHECK(testing::maxAbsDiff(robot.tipPosition(), tipBefore) == 0.0);
    CHECK(robot.configuration() == testing::kHomeConfig);
    CHECK(copy.configuration() == q2);
    CHECK(testing::maxAbsDiff(copy.tipPosition(), tipBefore) > 1e-6);
}

TEST_CASE("A moved-from robot self-heals on next use", "[unit][ctr]")
{
    CTR robot = testing::makeReferenceRobot();
    CTR other(std::move(robot));

    bvp_type guess{};
    REQUIRE(other.actuate(testing::kHomeConfig, guess));

    // NOLINTBEGIN(bugprone-use-after-move) — deliberately exercising the moved-from state
    bvp_type guess2{};
    const FKResult fk = robot.actuate(testing::kHomeConfig, guess2);
    REQUIRE(fk); // ensureSolver() must have recreated the solver — no crash, correct result
    CHECK(testing::maxAbsDiff(robot.tipPosition(), other.tipPosition()) < 1e-12);
    // NOLINTEND(bugprone-use-after-move)
}

TEST_CASE("Shape accessors are in meters and consistent", "[unit][ctr]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));

    const auto backbone = robot.backboneShape();
    REQUIRE_FALSE(backbone.empty());

    // Starts at the origin, ends at the tip, in METERS (no hidden ×1e3).
    CHECK(blaze::norm(backbone.front()) < 1e-12);
    CHECK(testing::maxAbsDiff(backbone.back(), robot.tipPosition()) < 1e-12);

    const auto shapes = robot.tubeShapes();
    REQUIRE_FALSE(shapes[0].empty());
    CHECK(testing::maxAbsDiff(shapes[0].back(), robot.tipPosition()) < 1e-12);
    // Inner tubes extend at least as far as outer ones (innermost-first ordering).
    CHECK(shapes[0].size() >= shapes[1].size());
    CHECK(shapes[1].size() >= shapes[2].size());

    // Tube 2/3 end near their distal arc-lengths: |r| <= arc length travelled.
    const auto distal = robot.distalEnds();
    CHECK(blaze::norm(shapes[1].back()) <= distal[1UL] + 1e-6);
    CHECK(blaze::norm(shapes[2].back()) <= distal[2UL] + 1e-6);
}

TEST_CASE("setTube swaps a tube and recomputes the segmentation", "[unit][ctr]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));
    const auto tipBefore = robot.tipPosition();

    // Stiffer innermost tube -> straighter robot -> different tip.
    TubeParams p = robot.tubes()[0UL].parameters();
    p.E *= 3.0;
    robot.setTube(0UL, Tube(p));

    bvp_type guess2{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess2));
    CHECK(testing::maxAbsDiff(robot.tipPosition(), tipBefore) > 1e-5);
}
