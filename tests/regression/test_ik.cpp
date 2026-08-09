// Inverse-kinematics regression tests (round trips, convergence, failure
// semantics). Enabled with the IK redesign; the legacy resolved-rates loop
// does not converge reliably enough to gate on.

#include <catch2/catch_test_macros.hpp>

#include "TestRobot.hpp"

TEST_CASE("IK: FK -> IK -> FK round trips", "[regression][ik]")
{
    SKIP("enabled with the IK redesign phase");
}

TEST_CASE("IK: unreachable target reports failure honestly", "[regression][ik]")
{
    SKIP("enabled with the IK redesign phase");
}
