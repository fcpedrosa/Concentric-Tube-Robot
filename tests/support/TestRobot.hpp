#pragma once

// Shared test fixtures: the canonical demo robot (matches examples/fk_ik_demo.cpp)
// and vector-comparison helpers.

#include "ctr/CTR.hpp"
#include "ctr/detail/mathOperations.hpp"

namespace ctr::testing
{

/// Young's modulus of Nitinol used by the reference robot [Pa].
inline constexpr double kE = 65.0e9;
/// Poisson ratio and derived shear modulus [Pa].
inline constexpr double kNu = 0.32;
inline constexpr double kG = kE / (2.0 * (1.0 + kNu));

/// The reference home configuration q = [β (m) | α (rad)].
inline const blaze::StaticVector<double, 6UL> kHomeConfig = {-120.0e-3, -100.0e-3, -80.0e-3, 0.0, 0.0, 0.0};

/// BVP tolerance used by all tests unless stated otherwise.
inline constexpr double kBVPTol = 1.0e-6;

/// The three reference tubes (innermost first) — identical to the demo.
[[nodiscard]] inline std::array<Tube, NUM_TUBES> makeReferenceTubes()
{
    return {
        Tube{{.OD = 0.92e-3, .ID = 0.80e-3, .E = kE, .G = kG, .ls = 190.0e-3, .lc = 60.0e-3,
              .u_ast = {1.0 / 40.0e-3, 0.0, 0.0}}},
        Tube{{.OD = 1.10e-3, .ID = 0.97e-3, .E = kE, .G = kG, .ls = 120.0e-3, .lc = 80.0e-3,
              .u_ast = {1.0 / 100.0e-3, 0.0, 0.0}}},
        Tube{{.OD = 1.40e-3, .ID = 1.20e-3, .E = kE, .G = kG, .ls = 90.0e-3, .lc = 40.0e-3,
              .u_ast = {1.0 / 140.0e-3, 0.0, 0.0}}},
    };
}

/// The canonical demo robot at its home configuration.
[[nodiscard]] inline CTR makeReferenceRobot(RootFindingMethod method = RootFindingMethod::ModifiedNewtonRaphson)
{
    return CTR(makeReferenceTubes(), kHomeConfig, kBVPTol, method);
}

/// A robot with the reference geometry but zero pre-curvature (perfectly straight tubes).
[[nodiscard]] inline CTR makeStraightRobot(RootFindingMethod method = RootFindingMethod::ModifiedNewtonRaphson)
{
    std::array<Tube, NUM_TUBES> tubes = {
        Tube{{.OD = 0.92e-3, .ID = 0.80e-3, .E = kE, .G = kG, .ls = 190.0e-3, .lc = 60.0e-3, .u_ast = {0.0, 0.0, 0.0}}},
        Tube{{.OD = 1.10e-3, .ID = 0.97e-3, .E = kE, .G = kG, .ls = 120.0e-3, .lc = 80.0e-3, .u_ast = {0.0, 0.0, 0.0}}},
        Tube{{.OD = 1.40e-3, .ID = 1.20e-3, .E = kE, .G = kG, .ls = 90.0e-3, .lc = 40.0e-3, .u_ast = {0.0, 0.0, 0.0}}},
    };
    return CTR(std::move(tubes), kHomeConfig, kBVPTol, method);
}

/// Max-norm distance between two 3-vectors.
[[nodiscard]] inline double maxAbsDiff(const blaze::StaticVector<double, 3UL> &a,
                                       const blaze::StaticVector<double, 3UL> &b)
{
    return blaze::linfNorm(a - b);
}

} // namespace ctr::testing
