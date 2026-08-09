// Cross-validation of the BVP solver suite and externally-loaded cases.
// Enabled once the solver fixes land (sanitizer no longer zeroes the moment
// unknowns; Broyden/dog-leg updates corrected).

#include <catch2/catch_test_macros.hpp>

#include "TestRobot.hpp"

TEST_CASE("All solvers agree on the converged shooting solution", "[regression][solvers]")
{
    SKIP("enabled with the solver-fixes phase");
}

TEST_CASE("Externally loaded robot converges (distal force)", "[regression][solvers][loaded]")
{
    SKIP("enabled with the solver-fixes phase");
}
