<div align="center">

### Concentric Tube Robot Kinematics — C++ Static Library

![CTR Robot](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/images/CTR_Assembly.png)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/LICENSE)
[![Build](https://github.com/fcpedrosa/Concentric-Tube-Robot/actions/workflows/ci.yml/badge.svg)](https://github.com/fcpedrosa/Concentric-Tube-Robot/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-doxygen-brightgreen.svg)](https://fcpedrosa.github.io/Concentric-Tube-Robot/html/index.html)

</div>

This repository contains a C++ static library for implementing the forward and inverse kinematics of a three‑tube concentric tube robot (CTR) based on a Cosserat theory kinematic model as described in the paper:

[1] D. C. Rucker, B. A. Jones, and R. J. Webster III, “A Geometrically Exact Model for Externally Loaded Concentric-Tube Continuum Robots,” IEEE Trans. Robot., vol. 26, no. 5, pp. 769–780, Oct. 2010, doi: 10.1109/TRO.2010.2062570.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Documentation](#documentation)
- [Usage](#usage)
- [GitHub Workflows](#github-workflows)
- [CMake Options](#cmake-options)
- [Folder Layout](#folder-layout)
- [Citing](#citing)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Introduction

This C++ static library provides a compact implementation of the forward and inverse kinematics of a three‑tube CTR. It targets research and prototyping workflows where reproducible results, modularity, and extensibility are needed.

## Features

- Forward kinematics for a three‑tube CTR.
- Inverse kinematics (differential position control).
- Modular tube/segment model with clear separation of concerns.
- Doxygen documentation with class diagrams and source browsing.

![Backbone Shape](https://github.com/fcpedrosa/Concentric-Tube-Robot/blob/main/images/Backbone.png)

## Installation

To use this library, you can either clone the repository and build it from source as a CMake project (recommended) or download the precompiled binaries for your platform from the [Releases](https://github.com/fcpedrosa/Concentric-Tube-Robot/releases) section.

## Building Requirements

Before using the library, ensure that you have the following libraries installed in your system:

* [Boost](https://www.boost.org/)
* [Blaze Library](https://bitbucket.org/blaze-lib/blaze/src/master/)
* [LAPACK](http://www.netlib.org/lapack/)
* [OpenBLAS](https://www.openblas.net/)
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

Local precompiled outputs are placed in the Releases folder whenever you build. After running a build, you should find the CTR library and the example executable in Releases. For GitHub-hosted binaries, use the [Releases](https://github.com/fcpedrosa/Concentric-Tube-Robot/releases) page.

## Documentation

You can generate documentation locally using the Doxygen target:

```bash
cmake --build build --target doc_doxygen
```

Open the output at `build/docs/html/index.html`.

## Usage

For detailed documentation on the library's API and usage, please refer to the [Documentation](https://fcpedrosa.github.io/Concentric-Tube-Robot/html/index.html) section.

## GitHub Workflows

This repository includes a GitHub Actions workflow that builds the library on every push and pull request. The build badge at the top of this README reflects the latest workflow status on the main branch.

To run the workflow manually:
1. Go to the Actions tab in your GitHub repository.
2. Select the build workflow.
3. Click Run workflow and choose the branch.

The workflow verifies that the project configures and builds successfully. If you want downloadable binaries, add an upload step to the workflow (for example, upload the build outputs from the build directory or the Releases folder) and attach them to a GitHub Release.

## CMake Options

- Uses C++23 by default.
- Builds a static library target `CTRlib` and an executable target `${PROJECT_NAME}`.
- Doxygen docs target: `doc_doxygen`.

## Folder Layout

- `ctr_library/include` — Public headers
- `ctr_library/src` — Library sources
- `executable` — Example executable
- `docs` — Generated documentation (if built)
- `images` — Figures used by docs and README

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