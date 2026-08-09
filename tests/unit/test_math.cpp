#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "ctr/detail/mathOperations.hpp"
#include <numbers>

using namespace ctr;
using Catch::Matchers::WithinAbs;

TEST_CASE("deg2Rad converts correctly", "[unit][math]")
{
    CHECK_THAT(math::deg2Rad(180.0), WithinAbs(std::numbers::pi, 1e-15));
    CHECK_THAT(math::deg2Rad(-90.0), WithinAbs(-std::numbers::pi / 2.0, 1e-15));
    STATIC_CHECK(math::deg2Rad(0.0) == 0.0);
}

TEST_CASE("euler2Quaternion + getSO3 reproduce a z-rotation", "[unit][math]")
{
    const double theta = 0.73;

    blaze::StaticVector<double, 4UL> h;
    math::euler2Quaternion(0.0, theta, 0.0, h); // attitude = rotation about z in this convention

    blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R;
    math::getSO3(h, R);

    const auto Rz = math::rotz(theta);
    for (std::size_t i = 0; i < 3UL; ++i)
        for (std::size_t j = 0; j < 3UL; ++j)
            CHECK_THAT(R(i, j), WithinAbs(Rz(i, j), 1e-12));
}

TEST_CASE("getSO3 accepts vector expressions and is self-normalizing", "[unit][math]")
{
    blaze::StaticVector<double, 8UL> padded{};
    blaze::StaticVector<double, 4UL> h;
    math::euler2Quaternion(0.0, 1.1, 0.0, h);
    h *= 3.7; // getSO3's 2/||h||^2 scaling must absorb the magnitude
    blaze::subvector<2UL, 4UL>(padded) = h;

    blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R;
    math::getSO3(blaze::subvector<2UL, 4UL>(padded), R); // view, no temporary

    const auto Rz = math::rotz(1.1);
    for (std::size_t i = 0; i < 3UL; ++i)
        for (std::size_t j = 0; j < 3UL; ++j)
            CHECK_THAT(R(i, j), WithinAbs(Rz(i, j), 1e-12));
}

TEST_CASE("hatOperator times vector equals the cross product", "[unit][math]")
{
    const blaze::StaticVector<double, 3UL> a{1.2, -0.7, 0.3};
    const blaze::StaticVector<double, 3UL> b{0.4, 2.1, -1.5};

    const blaze::StaticVector<double, 3UL> viaHat = math::hatOperator(a) * b;
    const blaze::StaticVector<double, 3UL> viaCross = blaze::cross(a, b);

    for (std::size_t i = 0; i < 3UL; ++i)
        CHECK_THAT(viaHat[i], WithinAbs(viaCross[i], 1e-14));
}

TEST_CASE("quaternionDiff matches 0.5*q*[0,u] quaternion product", "[unit][math]")
{
    // For unit quaternion h and curvature u (body frame): h' = 0.5 * h ⊗ [0, u].
    blaze::StaticVector<double, 4UL> h;
    math::euler2Quaternion(0.3, 0.9, -0.4, h);
    const blaze::StaticVector<double, 3UL> u{2.0, -1.0, 0.5};

    const auto hs = math::quaternionDiff(u, h);

    // Explicit quaternion product h ⊗ (0, u):
    const double w = h[0], x = h[1], y = h[2], z = h[3];
    const blaze::StaticVector<double, 4UL> expected{
        0.5 * (-x * u[0] - y * u[1] - z * u[2]), 0.5 * (w * u[0] + y * u[2] - z * u[1]),
        0.5 * (w * u[1] - x * u[2] + z * u[0]), 0.5 * (w * u[2] + x * u[1] - y * u[0])};

    for (std::size_t i = 0; i < 4UL; ++i)
        CHECK_THAT(hs[i], WithinAbs(expected[i], 1e-14));
}
