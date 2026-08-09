#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "ctr/Tube.hpp"
#include <numbers>
#include <stdexcept>

using namespace ctr;
using Catch::Matchers::WithinRel;

namespace
{
TubeParams referenceParams()
{
    return {.OD = 0.92e-3, .ID = 0.80e-3, .E = 65.0e9, .G = 24.6e9, .ls = 190.0e-3, .lc = 60.0e-3,
            .u_ast = {25.0, 0.0, 0.0}};
}
} // namespace

TEST_CASE("Tube stiffness follows the thin-tube formulas", "[unit][tube]")
{
    const TubeParams p = referenceParams();
    const Tube tube(p);

    const double od4 = p.OD * p.OD * p.OD * p.OD;
    const double id4 = p.ID * p.ID * p.ID * p.ID;
    const double I = std::numbers::pi / 64.0 * (od4 - id4);
    const double J = std::numbers::pi / 32.0 * (od4 - id4);

    CHECK_THAT(tube.getK(Stiffness::Bending), WithinRel(p.E * I, 1e-12));
    CHECK_THAT(tube.getK(Stiffness::Torsion), WithinRel(p.G * J, 1e-12));
}

TEST_CASE("Tube getters report the construction parameters", "[unit][tube]")
{
    const TubeParams p = referenceParams();
    const Tube tube(p);

    CHECK(tube.getTubeLength() == p.ls + p.lc);
    CHECK(tube.getStraightLen() == p.ls);
    CHECK(tube.getCurvLen() == p.lc);
    CHECK(tube.get_u_ast(0UL) == p.u_ast[0UL]);
    CHECK(tube.get_u_ast(1UL) == p.u_ast[1UL]);

    const TubeParams back = tube.parameters();
    CHECK(back.OD == p.OD);
    CHECK(back.ID == p.ID);
    CHECK(back.E == p.E);
    CHECK(back.G == p.G);
    CHECK(back.ls == p.ls);
    CHECK(back.lc == p.lc);
    CHECK(back.u_ast == p.u_ast);
}

TEST_CASE("Tube construction validates its parameters", "[unit][tube]")
{
    TubeParams p = referenceParams();

    SECTION("OD must exceed ID")
    {
        p.OD = p.ID;
        CHECK_THROWS_AS(Tube(p), std::invalid_argument);
    }
    SECTION("ID must be positive")
    {
        p.ID = 0.0;
        CHECK_THROWS_AS(Tube(p), std::invalid_argument);
    }
    SECTION("moduli must be positive")
    {
        p.E = 0.0;
        CHECK_THROWS_AS(Tube(p), std::invalid_argument);
    }
    SECTION("shear modulus must be positive")
    {
        p.G = -1.0;
        CHECK_THROWS_AS(Tube(p), std::invalid_argument);
    }
    SECTION("total length must be positive")
    {
        p.ls = 0.0;
        p.lc = 0.0;
        CHECK_THROWS_AS(Tube(p), std::invalid_argument);
    }
}
