// Analytic anchor tests (T1-T4): physics identities that hold independently of
// any golden values. These are the primary guards against kinematic
// regressions (sign errors, index swaps, frame mistakes).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestRobot.hpp"

#include <cmath>
#include <vector>

using namespace ctr;
using Catch::Matchers::WithinAbs;

// ─── T1: straight tubes ───────────────────────────────────────────────────────

TEST_CASE("T1: a robot of straight tubes is a straight line along z", "[regression][anchor]")
{
    CTR robot = testing::makeStraightRobot();

    auto checkStraight = [&](const blaze::StaticVector<double, 6UL> &q)
    {
        bvp_type guess{};
        REQUIRE(robot.actuate(q, guess));

        const double L1 = robot.tubes()[0UL].getTubeLength() + q[0UL]; // dist_end of tube 1
        const auto tip = robot.tipPosition();
        CHECK_THAT(tip[0UL], WithinAbs(0.0, 1e-9));
        CHECK_THAT(tip[1UL], WithinAbs(0.0, 1e-9));
        CHECK_THAT(tip[2UL], WithinAbs(L1, 1e-9));
    };

    SECTION("home configuration")
    {
        checkStraight(testing::kHomeConfig);
    }
    SECTION("arbitrary tube rotations do not bend a straight robot")
    {
        blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
        q[3UL] = 0.5;
        q[4UL] = -1.2;
        q[5UL] = 2.0;
        checkStraight(q);
    }
    SECTION("different retraction changes only the length")
    {
        blaze::StaticVector<double, 6UL> q = {-150.0e-3, -110.0e-3, -90.0e-3, 0.0, 0.0, 0.0};
        checkStraight(q);
    }
}

// ─── T2: torsionless planar robot vs closed-form arc chain ────────────────────

TEST_CASE("T2: equal-rotation planar robot matches the analytic arc chain", "[regression][anchor]")
{
    CTR robot = testing::makeReferenceRobot();
    bvp_type guess{};
    REQUIRE(robot.actuate(testing::kHomeConfig, guess));

    // With all alpha equal and pure-x pre-curvatures the robot is planar and
    // torsionless: each segment is a circular arc in the y-z plane with
    //   kappa_j = sum_curved(EI_i * ux_i) / sum_present(EI_i).
    const auto &tubes = robot.tubes();
    const blaze::StaticVector<double, NUM_TUBES> beta = robot.beta();

    std::vector<double> S{0.0};
    blaze::StaticVector<double, NUM_TUBES> distEnd, lenCurv, EI, ux;
    for (std::size_t i = 0; i < NUM_TUBES; ++i)
    {
        distEnd[i] = tubes[i].getTubeLength() + beta[i];
        lenCurv[i] = distEnd[i] - tubes[i].getCurvLen();
        EI[i] = tubes[i].getK(Stiffness::Bending);
        ux[i] = tubes[i].get_u_ast(0UL);
        S.push_back(lenCurv[i]);
        S.push_back(distEnd[i]);
    }
    std::sort(S.begin(), S.end());
    S.erase(std::unique(S.begin(), S.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), S.end());

    double phi = 0.0, ry = 0.0, rz = 0.0;
    for (std::size_t j = 0; j + 1UL < S.size(); ++j)
    {
        const double s0 = S[j], s1 = S[j + 1UL];
        const double mid = 0.5 * (s0 + s1);
        const double len = s1 - s0;

        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < NUM_TUBES; ++i)
        {
            if (distEnd[i] > mid) // tube present in this segment
            {
                den += EI[i];
                if (lenCurv[i] < mid) // segment lies in the tube's curved section
                    num += EI[i] * ux[i];
            }
        }
        const double kappa = num / den;

        if (std::fabs(kappa) < 1e-12)
        {
            ry += -std::sin(phi) * len;
            rz += std::cos(phi) * len;
        }
        else
        {
            ry += (std::cos(phi + kappa * len) - std::cos(phi)) / kappa;
            rz += (std::sin(phi + kappa * len) - std::sin(phi)) / kappa;
            phi += kappa * len;
        }
    }

    const auto tip = robot.tipPosition();
    CHECK_THAT(tip[0UL], WithinAbs(0.0, 1e-8));
    CHECK_THAT(tip[1UL], WithinAbs(ry, 2e-5));
    CHECK_THAT(tip[2UL], WithinAbs(rz, 2e-5));
}

// ─── T3: twist-rate identity θ'ᵢ = u_zi − u_z1 ───────────────────────────────

TEST_CASE("T3: recorded twist angles satisfy theta' = uz_i - uz_1", "[regression][anchor]")
{
    using namespace StateIdx;

    CTR robot = testing::makeReferenceRobot();
    blaze::StaticVector<double, 6UL> q = testing::kHomeConfig;
    q[4UL] = math::deg2Rad(60.0);
    q[5UL] = math::deg2Rad(120.0);

    bvp_type guess{};
    REQUIRE(robot.actuate(q, guess));

    const auto states = robot.states();
    const auto s = robot.arcLengthSamples();
    REQUIRE(states.size() == s.size());
    REQUIRE(states.size() > 10UL);

    const auto distal = robot.distalEnds();

    std::size_t checked = 0UL;
    for (std::size_t k = 0; k + 1UL < states.size(); ++k)
    {
        const double ds = s[k + 1UL] - s[k];
        if (ds < 1e-6) // duplicate sample at a segment boundary
            continue;

        for (std::size_t tube = 1UL; tube < NUM_TUBES; ++tube)
        {
            // The identity θ'ᵢ = uzᵢ − uz₁ holds only while tube i exists;
            // beyond its distal end the ODE freezes θᵢ and uzᵢ by construction.
            if (s[k + 1UL] > distal[tube] - 1e-6)
                continue;

            const double thetaRate = (states[k + 1UL][THETA_1 + tube] - states[k][THETA_1 + tube]) / ds;
            const double uzDiff = 0.5 * ((states[k][UZ_1 + tube] - states[k][UZ_1]) +
                                         (states[k + 1UL][UZ_1 + tube] - states[k + 1UL][UZ_1]));

            // Trapezoidal comparison: theta' must equal uz_i - uz_1, NOT theta_i - theta_1.
            CHECK_THAT(thetaRate, WithinAbs(uzDiff, 1e-3));
            ++checked;
        }
    }
    // The identity must actually have been exercised (non-planar configuration).
    CHECK(checked > 20UL);
}

// ─── T4: rotational equivariance ─────────────────────────────────────────────

TEST_CASE("T4: rotating all tubes together rotates the tip rigidly", "[regression][anchor]")
{
    CTR robot = testing::makeReferenceRobot();

    blaze::StaticVector<double, 6UL> qa = testing::kHomeConfig;
    qa[3UL] = 0.3;
    qa[4UL] = 1.0;
    qa[5UL] = -0.6;

    bvp_type guessA{};
    REQUIRE(robot.actuate(qa, guessA));
    const auto tipA = robot.tipPosition();

    const double delta = 0.7;
    blaze::StaticVector<double, 6UL> qb = qa;
    qb[3UL] += delta;
    qb[4UL] += delta;
    qb[5UL] += delta;

    bvp_type guessB{};
    REQUIRE(robot.actuate(qb, guessB));
    const auto tipB = robot.tipPosition();

    const blaze::StaticVector<double, 3UL> expected = math::rotz(delta) * tipA;
    CHECK(testing::maxAbsDiff(tipB, expected) < 1e-7);
}
