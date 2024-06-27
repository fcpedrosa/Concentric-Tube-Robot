<div align="center">

### Concentric Tube Robot Kinematics -- C++ Static Library

![CTR Robot](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/images/CTR_Assembly.png)

</div>


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/LICENSE)
[![Build Status](https://travis-ci.com/yourusername/repo-name.svg?branch=main)](https://github.com/fcpedrosa/Concentric-Tube-Robot)

This repository contains a C++ static library for implementing the forward and inverse kinematics of a three-tube concentric tube robot (CTR) based on the Cosserat Theory-based kinematic model as described in the paper:

[1] D. C. Rucker, B. A. Jones, and R. J. Webster III, “A Geometrically Exact Model for Externally Loaded Concentric-Tube Continuum Robots,” IEEE Trans. Robot., vol. 26, no. 5, pp. 769–780, Oct. 2010, doi: 10.1109/TRO.2010.2062570.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Documentation & Usage](#usage)
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
- Comprehensive [documentation](https://fcpedrosa.github.io/Concentric-Tube-Robot/html/index.html) explaining the kinematic model and usage of the library.

![Backbone Shape](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/images/Backbone.png)

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

Download the precompiled binaries for your platform from the [Releases](https://github.com/fcpedrosa/Concentric-Tube-Robot/releases) section (TBA). Add the library to your C++ project's dependencies and include the appropriate header files.

## Usage

For detailed documentation on the library's API and usage, please refer to the [Documentation](https://fcpedrosa.github.io/Concentric-Tube-Robot/html/index.html) section.

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