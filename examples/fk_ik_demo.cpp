// Demonstrates the intended API usage of the CTR kinematics library:
// construct three tubes, assemble a robot, run forward kinematics (FK),
// then drive the tip to a Cartesian target with inverse kinematics (IK).

#include <chrono>
#include <iostream>
#include "CTR.hpp"

using namespace ctr;

int main()
{
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

    // -- ** -- Precurvature vectors (for curved portions of the tubes) -- ** -- [u_x* u_y* 0]
    constexpr blaze::StaticVector<double, 3UL> u1 = {1.00 / R1, 0.00, 0.00};
    constexpr blaze::StaticVector<double, 3UL> u2 = {1.00 / R2, 0.00, 0.00};
    constexpr blaze::StaticVector<double, 3UL> u3 = {1.00 / R3, 0.00, 0.00};

    // --** --Lengths of the tubes' straight sections [m] -- ** --
    constexpr blaze::StaticVector<double, 3UL> ls = {190.00E-3, 120.00E-3, 90.00E-3};

    // --** --Lengths of the tubes' curved sections [m] -- ** --
    constexpr blaze::StaticVector<double, 3UL> lc = {60.00E-3, 80.00E-3, 40.00E-3};

    // --** --Outer and Inner diameters of the tubes [m]--** --
    constexpr blaze::StaticVector<double, 3UL> OD = {0.92e-3, 1.10E-3, 1.40E-3};
    constexpr blaze::StaticVector<double, 3UL> ID = {0.80E-3, 0.97e-3, 1.20E-3};

    // # # # # # ---- Instantiating the three Tube objects ---- # # # # #
    std::shared_ptr<Tube> T1 = std::make_shared<Tube>(OD[0UL], ID[0UL], E, G, ls[0UL], lc[0UL], u1); // innermost tube
    std::shared_ptr<Tube> T2 =
        std::make_shared<Tube>(OD[1UL], ID[1UL], E, G, ls[1UL], lc[1UL], u2); // intermediate tube
    std::shared_ptr<Tube> T3 = std::make_shared<Tube>(OD[2UL], ID[2UL], E, G, ls[2UL], lc[2UL], u3); // outermost tube

    // instantiating an array of smart pointers to CTR component tubes
    std::array<std::shared_ptr<Tube>, 3UL> Tb = {T1, T2, T3};

    // initial joint actuation values "home position" - q = [Beta Alpha]
    blaze::StaticVector<double, 3UL> Beta_0 = {-120.00E-3, -100.00E-3, -80.00E-3};
    blaze::StaticVector<double, 3UL> Alpha_0 = {mathOp::deg2Rad(0.00), mathOp::deg2Rad(0.00), mathOp::deg2Rad(0.00)};

    blaze::StaticVector<double, 6UL> q_0;
    blaze::subvector<0UL, 3UL>(q_0) = Beta_0;
    blaze::subvector<3UL, 3UL>(q_0) = Alpha_0;

    // Determining the accuracy of BVP solutions
    constexpr double Tol = 1.00E-6;

    // tolerance for position control (0.5 mm)
    constexpr double pos_tol = 5.00E-4;

    // # # # # # ---- Instantiating the CTR object ---- # # # # #
    CTR CTR_robot(Tb, q_0, Tol, mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON);

    // initial guess for the BVP
    bvp_type initGuess;

    // ************************** Actuating the CTR and solving the corresponding BVP **************************
    auto start = std::chrono::high_resolution_clock::now();
    const bool fk_ok = CTR_robot.actuate_CTR(initGuess, q_0);
    auto finish = std::chrono::high_resolution_clock::now();

    auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "FK " << (fk_ok ? "converged" : "FAILED to converge") << " in " << elapsed_us
              << " microseconds. Tip position [m]:\n"
              << CTR_robot.getTipPos() << std::endl;

    // >>> Driving the CTR tip to a Cartesian target
    blaze::StaticVector<double, 3UL> target = {-0.053210, 0.043606, 0.179527}, tip_pos;

    // inverse kinematics
    start = std::chrono::high_resolution_clock::now();
    const bool ik_ok = CTR_robot.posCTRL(initGuess, target, pos_tol);
    finish = std::chrono::high_resolution_clock::now();

    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "IK " << (ik_ok ? "succeeded" : "FAILED") << " in " << elapsed_ms << " milliseconds.\n" << std::endl;

    tip_pos = CTR_robot.getTipPos();

    std::cout << "Target [m]: " << blaze::trans(target) << "Tip position [m]: " << blaze::trans(tip_pos)
              << "Joint values (IK solution): " << blaze::trans(CTR_robot.getConfiguration())
              << "Position error: " << blaze::norm(tip_pos - target) * 1.00E3 << " [mm]" << std::endl;

    return (fk_ok && ik_ok) ? 0 : 1;
}
