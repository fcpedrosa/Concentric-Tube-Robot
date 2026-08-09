// Golden-value FK regression tests.
//
// Golden values are captured AFTER the integration-strategy phase (fixed-step
// RK4) so they are not invalidated by planned numerical changes. Until then
// these tests are skipped.

#include <catch2/catch_test_macros.hpp>

#include "TestRobot.hpp"

TEST_CASE("Golden FK: 20-configuration grid", "[regression][golden]")
{
    SKIP("golden values are captured in the integration-strategy phase");
}

TEST_CASE("Golden FK: loaded configurations", "[regression][golden]")
{
    SKIP("golden values are captured in the integration-strategy phase");
}
