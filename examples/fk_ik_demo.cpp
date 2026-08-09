// Demonstrates the intended API usage of the CTR kinematics library:
// construct three tubes, assemble a robot, run forward kinematics (FK),
// then drive the tip to a Cartesian target with inverse kinematics (IK).

#include <chrono>
#include <iostream>
#include "ctr/CTR.hpp"
#include "ctr/detail/mathOperations.hpp"

int main()
{
    using namespace ctr;

    //  # # # # # # # # ---- Properties of Nitinol Tubes ---- # # # # # # # #
    // Young's modulus [Pa]
    constexpr double E = 65.00E9;
    // Poisson's ratio
    constexpr double nu = 0.32;
    // Shear modulus [Pa]
    constexpr double G = E / (2.00 * (1.00 + nu));

    // Precurvature radii for the tubes [m]
    constexpr double R1 = 40.000E-3; // (4 cm curvature radius)
    constexpr double R2 = 100.00E-3; // (10 cm curvature radius)
    constexpr double R3 = 140.00E-3; // (14 cm curvature radius)

    // # # # # # ---- The three component tubes, innermost first ---- # # # # #
    // All quantities SI: meters, Pascals, 1/m.
    std::array<Tube, NUM_TUBES> tubes = {
        Tube{{.OD = 0.92E-3, .ID = 0.80E-3, .E = E, .G = G, .ls = 190.00E-3, .lc = 60.00E-3,
              .u_ast = {1.00 / R1, 0.00, 0.00}}}, // innermost tube
        Tube{{.OD = 1.10E-3, .ID = 0.97E-3, .E = E, .G = G, .ls = 120.00E-3, .lc = 80.00E-3,
              .u_ast = {1.00 / R2, 0.00, 0.00}}}, // intermediate tube
        Tube{{.OD = 1.40E-3, .ID = 1.20E-3, .E = E, .G = G, .ls = 90.00E-3, .lc = 40.00E-3,
              .u_ast = {1.00 / R3, 0.00, 0.00}}}, // outermost tube
    };

    // initial joint actuation values "home position" - q = [beta | alpha]
    const blaze::StaticVector<double, 6UL> q_0 = {-120.00E-3, -100.00E-3, -80.00E-3, // beta [m]
                                                  0.00, 0.00, 0.00};                 // alpha [rad]

    // BVP tolerance and IK position tolerance (0.5 mm)
    constexpr double Tol = 1.00E-6;
    constexpr double pos_tol = 5.00E-4;

    // # # # # # ---- Instantiating the CTR object ---- # # # # #
    CTR robot(tubes, q_0, Tol, RootFindingMethod::ModifiedNewtonRaphson);

    // initial guess for the BVP (zero: cold start)
    bvp_type initGuess;

    // ************************ Forward kinematics ************************
    auto start = std::chrono::high_resolution_clock::now();
    const FKResult fk = robot.actuate(q_0, initGuess);
    auto finish = std::chrono::high_resolution_clock::now();

    const auto fk_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "FK " << (fk ? "converged" : "FAILED to converge") << " in " << fk_us << " us ("
              << fk.iterations << " iterations, residual " << fk.residual << ").\n"
              << "Tip position [m]: " << blaze::trans(robot.tipPosition()) << std::endl;

    // ************************ Inverse kinematics ************************
    // Build a verified-reachable target: run FK at a reference configuration,
    // take its tip position, then return to the home configuration and ask IK
    // to steer the tip back to that target.
    // (The historical demo target {-0.0532, 0.0436, 0.1795} lies OUTSIDE this
    // tube set's reachable workspace — multi-start IK bottoms out ~26 mm away.
    // The legacy solver claimed success on it because its return value
    // reported BVP convergence, not target attainment.)
    blaze::StaticVector<double, 6UL> q_ref = q_0;
    q_ref[4UL] = math::deg2Rad(60.0);
    q_ref[5UL] = math::deg2Rad(-45.0);
    std::ignore = robot.actuate(q_ref, initGuess);
    const blaze::StaticVector<double, 3UL> target = robot.tipPosition();
    std::ignore = robot.actuate(q_0, initGuess); // back to home

    start = std::chrono::high_resolution_clock::now();
    const IKResult ik = robot.solveIK(target, pos_tol, initGuess);
    finish = std::chrono::high_resolution_clock::now();

    const auto ik_ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "IK " << (ik ? "reached the target" : "FAILED to reach the target") << " in " << ik_ms << " ms ("
              << ik.iterations << " iterations).\n"
              << std::endl;

    std::cout << "Target [m]:       " << blaze::trans(target) << "Tip position [m]: "
              << blaze::trans(robot.tipPosition()) << "Joint values (IK solution): " << blaze::trans(ik.q)
              << "Position error: " << ik.positionError * 1.00E3 << " [mm]" << std::endl;

    return (fk && ik) ? 0 : 1;
}
