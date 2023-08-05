# Concentric Tube Robot Kinematics -- C++ Static Library

![CTR Robot](https://drive.google.com/file/d/1JG_TEgWo15-uIoW0uRBfLdJjFkRGMBR6/view?usp=drive_link)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/LICENSE)
[![Build Status](https://travis-ci.com/yourusername/repo-name.svg?branch=main)](https://github.com/fcpedrosa/Concentric-Tube-Robot)

This repository contains a C++ static library for implementing the forward and inverse kinematics of a three-tube concentric tube robot (CTR) based on the Cosserat Theory-based kinematic model as described in the paper:

[1] D. C. Rucker, B. A. Jones, and R. J. Webster III, “A Geometrically Exact Model for Externally Loaded Concentric-Tube Continuum Robots,” IEEE Trans. Robot., vol. 26, no. 5, pp. 769–780, Oct. 2010, doi: 10.1109/TRO.2010.2062570.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Introduction

This C++ static library provides a barebone implementation of the forward and inverse kinematics of a three-tube concentric tube robot based on the Cosserat Theory-based kinematic model. As part of my PhD at Western University in London, Ontario, Canada, I have developed this C++ static library to provide a reliable and efficient implementation of the forward and inverse kinematics for a three-tube CTR. The library is based on the Cosserat Theory-based kinematic model as described in [1]. Concentric tube robots have applications in medical robotics and minimally invasive procedures.

## Features

- Implementation of the forward kinematics of a three-tube CTR.
- Implementation of the inverse kinematics of a three-tube CTR.
- Easy-to-use API for integrating the library into your C++ projects.
- Comprehensive documentation explaining the kinematic model and usage of the library.

## Installation

To use this library, you can either clone the repository and build it from source as a CMAKE project (recommended) or download the precompiled binaries for your platform from the [Releases](https://github.com/fcpedrosa/Concentric-Tube-Robot/releases) section.

## Building Requirements

Before using the library, ensure that you have the following libraries installed in your system:

* [Boost](https://www.boost.org/)
* [Blaze Library](https://bitbucket.org/blaze-lib/blaze/src/master/)
* [LAPACK](http://www.netlib.org/lapack/)
* [openBLAS](https://www.openblas.net/)
* [TBB (Threading Building Blocks)](https://www.threadingbuildingblocks.org/)

### Building from Source

1. Clone the repository:

```bash
git clone https://github.com/fcpedrosa/Concentric-Tube-Robot
```

2. Build the library:

```bash
cd repo-name
mkdir build && cd build
cmake ..
make -j4
```

3. The library will be built as a static library (`CTR.a`).

### Using Precompiled Binaries

Download the precompiled binaries for your platform from the [Releases](https://github.com/fcpedrosa/Concentric-Tube-Robot/releases) section. Add the library to your C++ project's dependencies and include the appropriate header files.

## Usage

To use the library in your C++ project, you need to include the necessary header files and link against the library. Here's a simple example:

```cpp
#include <iostream>
// include the CTR library header
#include "CTR.hpp"

int main()
{	
	//  # # # # # # # # ---- Properties of Nitinol Tubes ---- # # # # # # # #
	double E = 58.00E9; // Young's modulus GPa
	double G = 25.50E9; // Shear modulus GPa

	// Precurvature radii for the tubes
	double R1 = 0.04; // (4cm curvature radius)
	double R2 = 0.10; // (10 cm curvature radius)
	double R3 = 0.14;  // (14 cm curvature radius)

	// -- ** -- Precurvature vectors (for curved portions of the tubes) -- ** -- [u_x* u_y* 0]
	blaze::StaticVector<double, 3UL> u1, u2, u3;
	u1 = {1.00 / R1, 0.00, 0.00};
	u2 = {1.00 / R2, 0.00, 0.00};
	u3 = {1.00 / R3, 0.00, 0.00};

    // defining the kinematic parameters of the robot
	blaze::StaticVector<double, 3UL> ls, lc, OD, ID;
	// --** --Lengths of the tubes' straight sections (meters) -- ** --
	ls = {190.00E-3, 120.00E-3, 90.00E-3}; // 190, 120, 100

	// --** --Lengths of the tubes' curved sections (meters) -- ** --
	lc = {60.00E-3, 80.00E-3, 40.00E-3}; // 60, 80, 50;

	// --** --Outer and Inner diameters of the tubes (meters)--** --
	OD = {0.92e-3, 1.10E-3, 1.40E-3};
	ID = {0.80E-3, 0.97e-3, 1.20E-3};

	// # # # # # ---- Instantiating the three Tube objects ---- # # # # #
	std::shared_ptr<Tube> T1 = std::make_shared<Tube>(OD[0UL], ID[0UL], E, G, ls[0UL], lc[0UL], u1); // innermost tube
	std::shared_ptr<Tube> T2 = std::make_shared<Tube>(OD[1UL], ID[1UL], E, G, ls[1UL], lc[1UL], u2); // intermediate tube
	std::shared_ptr<Tube> T3 = std::make_shared<Tube>(OD[2UL], ID[2UL], E, G, ls[2UL], lc[2UL], u3); // outermost tube

	// instantiating an array of smart pointers to CTR component tubes
	std::array<std::shared_ptr<Tube>, 3UL> Tb = {T1, T2, T3};

	// initial joint actuation values "home position" - q = [Beta Alpha]
	blaze::StaticVector<double, 3UL> Beta_0 = {-120.00E-3, -100.00E-3, -80.00E-3}; // 130, 100, 50 | 130, 100, 
	blaze::StaticVector<double, 3UL> Alpha_0 = {mathOp::deg2Rad(0.00), mathOp::deg2Rad(0.00), mathOp::deg2Rad(0.00)};

	blaze::StaticVector<double, 6UL> q_0;
	blaze::subvector<0UL, 3UL>(q_0) = Beta_0;
	blaze::subvector<3UL, 3UL>(q_0) = Alpha_0;

	// Determining the accuracy of BVP solutions
	double Tol = 1.00E-6;

	// tolerance for position control (0.5 mm)
	double pos_tol = 5.00E-4;

	// # # # # # ---- Instantiating the CTR object ---- # # # # #
	CTR CTR_robot(Tb, q_0, Tol, mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON);

	// initial guess for the BVP
	blaze::StaticVector<double, 5UL> initGuess;

	// Actuates the robot to the configuration q_0 and solves the corresponding FK problem
	CTR_robot.actuate_CTR(initGuess, q_0);	

	// Testing the differential kinematics-based position control for the CTR
	blaze::StaticVector<double, 3UL> target = { -0.053210, 0.043606, 0.179527 }, tip_pos;

	// inverse kinematics
	CTR_robot.posCTRL(initGuess, target, pos_tol);
	
	tip_pos = CTR_robot.getTipPos();
	
	std::cout << "Target is: " << blaze::trans(target) 
			  << "Tip position is: " << blaze::trans(tip_pos) 
			  << "Joint values (IK solution): " << blaze::trans(CTR_robot.getConfiguration()) 
			  << "Position error: " << blaze::norm( tip_pos - target ) 
			  << std::endl;
	
	return 0;
}
```

For detailed documentation on the library's API and usage, refer to the [Documentation](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/docs/README.md) section.

## Examples

The [examples](https://github.com/fcpedrosa/Concentric-Tube-Robot/tree/main/examples) directory contains some usage examples to help you get started.

## Contributing

Contributions to this project are welcome! If you find any issues or have ideas for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/LICENSE) file for details.

## Contact

For any questions or inquiries, feel free to contact me:

Filipe Pedrosa
Email: fpedrosa@uwo.ca

Please note that this project is a research-oriented implementation and comes with no warranty or support. Use it at your own risk.