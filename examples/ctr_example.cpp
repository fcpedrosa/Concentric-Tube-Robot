#include <array>
#include <iostream>
#include <memory>

#include <blaze/Math.h>

#include "CTR.hpp"
#include "Tube.hpp"
#include "mathOperations.hpp"

int main()
{
    // Properties of Nitinol Tubes
    constexpr double E = 58.00E9; // Young's modulus [Pa]
    constexpr double nu = 0.32;   // Poisson's ratio
    constexpr double G = E / (2.00 * (1.00 + nu)); // Shear modulus [Pa]

    // Precurvature radii for the tubes
    constexpr double R1 = 40.00E-3;  // 4 cm curvature radius
    constexpr double R2 = 100.00E-3; // 10 cm curvature radius
    constexpr double R3 = 140.00E-3; // 14 cm curvature radius

    // Precurvature vectors for the curved portions of the tubes
    constexpr blaze::StaticVector<double, 3UL> u1 = {1.00 / R1, 0.00, 0.00};
    constexpr blaze::StaticVector<double, 3UL> u2 = {1.00 / R2, 0.00, 0.00};
    constexpr blaze::StaticVector<double, 3UL> u3 = {1.00 / R3, 0.00, 0.00};

    // Kinematic parameters of the robot
    constexpr blaze::StaticVector<double, 3UL> ls = {190.00E-3, 120.00E-3, 90.00E-3}; // straight lengths
    constexpr blaze::StaticVector<double, 3UL> lc = {60.00E-3, 80.00E-3, 40.00E-3};   // curved lengths
    constexpr blaze::StaticVector<double, 3UL> OD = {0.92e-3, 1.10E-3, 1.40E-3};      // outer diameters
    constexpr blaze::StaticVector<double, 3UL> ID = {0.80E-3, 0.97e-3, 1.20E-3};      // inner diameters

    // Instantiate Tube objects
    auto T1 = std::make_shared<Tube>(OD[0UL], ID[0UL], E, G, ls[0UL], lc[0UL], u1); // innermost tube
    auto T2 = std::make_shared<Tube>(OD[1UL], ID[1UL], E, G, ls[1UL], lc[1UL], u2); // middle tube
    auto T3 = std::make_shared<Tube>(OD[2UL], ID[2UL], E, G, ls[2UL], lc[2UL], u3); // outermost tube

    std::array<std::shared_ptr<Tube>, 3UL> tubes = {T1, T2, T3};

    // Initial joint actuation values "home position" - q = [Beta Alpha]
    blaze::StaticVector<double, 3UL> Beta_0 = {-120.00E-3, -100.00E-3, -80.00E-3};
    blaze::StaticVector<double, 3UL> Alpha_0 = {
        mathOp::deg2Rad(0.00),
        mathOp::deg2Rad(0.00),
        mathOp::deg2Rad(0.00)};

    blaze::StaticVector<double, 6UL> q_0;
    blaze::subvector<0UL, 3UL>(q_0) = Beta_0;
    blaze::subvector<3UL, 3UL>(q_0) = Alpha_0;

    // BVP accuracy and position control tolerance
    const double bvp_tol = 1.00E-6;
    const double pos_tol = 5.00E-4; // 0.5 mm

    // Instantiate CTR
    CTR ctr(tubes, q_0, bvp_tol, mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON);

    // Initial guess for the BVP
    blaze::StaticVector<double, 5UL> initGuess{};

    // External wrench at distal end (optional)
    blaze::StaticVector<double, 3UL> moment = {0.00, -0.05, 0.00};
    blaze::StaticVector<double, 3UL> force = {0.00, -0.15, -0.55};
    ctr.setDistalMoment(moment);
    ctr.setDistalForce(force);

    // Forward kinematics: actuate the robot to q_0
    if (!ctr.actuate_CTR(initGuess, q_0))
    {
        std::cerr << "FK did not converge." << std::endl;
        return 1;
    }

    std::cout << "Tip position (FK): " << blaze::trans(ctr.getTipPos()) << std::endl;

    // Inverse kinematics: steer to a target point
    blaze::StaticVector<double, 3UL> target = {-0.053210, 0.043606, 0.179527};
    if (!ctr.posCTRL(initGuess, target, pos_tol))
    {
        std::cerr << "IK did not converge." << std::endl;
        return 1;
    }

    auto tip_pos = ctr.getTipPos();
    std::cout << "Target: " << blaze::trans(target)
              << "Tip: " << blaze::trans(tip_pos)
              << "q (IK): " << blaze::trans(ctr.getConfiguration())
              << "Error [mm]: " << blaze::norm(tip_pos - target) * 1.00E3
              << std::endl;

    return 0;
}
