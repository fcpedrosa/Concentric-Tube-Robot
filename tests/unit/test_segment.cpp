#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestRobot.hpp"
#include "ctr/Segment.hpp"

#include <algorithm>

using namespace ctr;
using Catch::Matchers::WithinAbs;

TEST_CASE("Segment decomposition of the reference robot at home", "[unit][segment]")
{
    const auto tubes = testing::makeReferenceTubes();
    const blaze::StaticVector<double, NUM_TUBES> beta = {-120.0e-3, -100.0e-3, -80.0e-3};

    const Segment seg(tubes, beta);
    const auto &S = seg.getTransitionPoints();

    // dist_end = L_i + beta_i = [0.13, 0.10, 0.05]; len_curv = dist_end - lc = [0.07, 0.02, 0.01]
    REQUIRE(S.size() == 7UL);
    CHECK(std::is_sorted(S.begin(), S.end()));
    CHECK_THAT(S.front(), WithinAbs(0.0, 1e-12));
    CHECK_THAT(S.back(), WithinAbs(0.13, 1e-9));

    const std::array<double, 7UL> expected = {0.0, 0.01, 0.02, 0.05, 0.07, 0.10, 0.13};
    for (std::size_t i = 0; i < 7UL; ++i)
        CHECK_THAT(S[i], WithinAbs(expected[i], 1e-9));

    const auto &dist = seg.getDistalEnds();
    CHECK_THAT(dist[0UL], WithinAbs(0.13, 1e-9));
    CHECK_THAT(dist[1UL], WithinAbs(0.10, 1e-9));
    CHECK_THAT(dist[2UL], WithinAbs(0.05, 1e-9));
}

TEST_CASE("Segment clamps transition points to non-negative arc length", "[unit][segment]")
{
    const auto tubes = testing::makeReferenceTubes();
    // Tube 1 retracted so far that its curved section starts inside the
    // actuation unit: len_curv_1 = 0.25 - 0.21 - 0.06 = -0.02 -> clamped to 0.
    const blaze::StaticVector<double, NUM_TUBES> beta = {-210.0e-3, -100.0e-3, -80.0e-3};

    const Segment seg(tubes, beta);
    const auto &S = seg.getTransitionPoints();

    REQUIRE_FALSE(S.empty());
    CHECK(S.front() == 0.0);
    for (const double s : S)
        CHECK(s >= 0.0);
    CHECK(std::is_sorted(S.begin(), S.end()));

    const auto &dist = seg.getDistalEnds();
    for (std::size_t i = 0; i < NUM_TUBES; ++i)
        CHECK(dist[i] >= 0.0);
}

TEST_CASE("Segment stiffness/pre-curvature matrices honor tube presence", "[unit][segment]")
{
    const auto tubes = testing::makeReferenceTubes();
    const blaze::StaticVector<double, NUM_TUBES> beta = {-120.0e-3, -100.0e-3, -80.0e-3};

    const Segment seg(tubes, beta);
    const auto &EI = seg.get_EI();
    const auto &U_x = seg.get_U_x();

    REQUIRE(EI.columns() == 6UL);

    // Segments: [0,.01] [.01,.02] [.02,.05] [.05,.07] [.07,.10] [.10,.13]
    // Tube 1 (dist_end 0.13) spans all segments.
    for (std::size_t j = 0; j < 6UL; ++j)
        CHECK(EI(0UL, j) > 0.0);

    // Tube 3 (dist_end 0.05) is present only in segments 0-2.
    for (std::size_t j = 0; j < 3UL; ++j)
        CHECK(EI(2UL, j) > 0.0);
    for (std::size_t j = 3UL; j < 6UL; ++j)
        CHECK(EI(2UL, j) == 0.0);

    // Tube 1 curved section starts at 0.07 -> pre-curvature only in segments 4-5.
    for (std::size_t j = 0; j < 4UL; ++j)
        CHECK(U_x(0UL, j) == 0.0);
    for (std::size_t j = 4UL; j < 6UL; ++j)
        CHECK(U_x(0UL, j) > 0.0);

    // Tube 3 curved section spans [0.01, 0.05] -> segments 1-2.
    CHECK(U_x(2UL, 0UL) == 0.0);
    CHECK(U_x(2UL, 1UL) > 0.0);
    CHECK(U_x(2UL, 2UL) > 0.0);
}
