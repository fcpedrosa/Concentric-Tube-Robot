# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C++23 static library implementing forward and inverse kinematics of a three-tube concentric tube robot (CTR), based on the Cosserat-theory kinematic model from Rucker, Jones & Webster, *IEEE Trans. Robot.* 2010. Research code developed for medical/continuum robotics.

## Build

Requires CMake ≥ 3.26, a C++23 compiler, and system-installed dependencies: Boost (odeint), Blaze, OpenBLAS, LAPACK, and TBB.

```bash
cmake -S . -B build            # defaults to Release if no build type given
cmake --build build -j
./build/bin/CTR_General_Use    # demo executable
```

Outputs: static library `build/lib/libCTRlib.a` (target `CTRlib`), executable `build/bin/CTR_General_Use` (target `ctr_general_use`).

Useful CMake options (all default OFF):
- `-DCTR_ENABLE_WARNINGS=ON` — -Wall -Wextra -Wpedantic (add `-DCTR_WARNINGS_AS_ERRORS=ON` for -Werror)
- `-DCTR_ENABLE_SANITIZERS=ON` — ASan + UBSan

There are no unit tests: `include(CTest)` is present but no `add_test()` targets exist. Verify changes by building and running the demo executable (`executable/Main.cpp`), which prints tip position, FK/IK timings, and IK position error.

## Architecture

Two targets: the library in `ctr_library/` (headers in `include/`, sources in `src/`) and a demo driver in `executable/Main.cpp` showing the intended API usage. All library code lives in namespace `ctr`, except free math utilities in namespace `mathOp` (`mathOperations.hpp`).

### The kinematics pipeline (the big picture)

Forward kinematics is a boundary value problem solved by a shooting method. The flow, orchestrated by the `CTR` class (`CTR.cpp`):

1. `CTR::actuate_CTR(initGuess, q)` sets the joint configuration and delegates to the active `BVPSolver` strategy.
2. The solver iteratively calls `CTR::ODESolver(initGuess)`: it resets state, asks `Segment` for per-segment parameters, then integrates the 15-dim ODE state segment-by-segment along the backbone arc-length with an 8th-order adaptive Adams-Bashforth-Moulton stepper (Boost.Odeint), returning the 5-dim residue of the distal boundary conditions.
3. The root-finder drives that residue below `m_accuracy` using the finite-difference `CTR::jac_BVP`.
4. Inverse kinematics (`CTR::posCTRL`) wraps FK in a resolved-rates loop using the finite-difference 3×6 kinematic Jacobian (`CTR::jacobian`).

### Key components

- **`CTRTypes.hpp`** — compile-time constants (`NUM_TUBES = 3`, `BVP_DIM = 5`), the `state_type` (15-dim) / `bvp_type` (5-dim) aliases, and `StateIdx` named indices into the ODE state vector. Use `StateIdx` constants, not raw indices.
- **`Tube`** — geometry and material of one tube. Derived stiffnesses (EI, GJ) are computed on demand by `getK()`, never cached. Selected via the `Stiffness` enum.
- **`Segment`** — splits the backbone into arc-length segments at tube transition points and caches per-segment stiffness/pre-curvature matrices. `recalculateSegments()` is private; `CTR` is a `friend` and is the only class allowed to trigger recomputation (must happen whenever linear actuation β changes).
- **`ODESystem`** — Boost.Odeint-compatible functor implementing the ODE right-hand side; reparameterized via `setEquationParameters()` before each segment.
- **`BVPSolver`** (`BVPSolver.hpp/.cpp`) — strategy interface with six concrete root-finders (Newton-Raphson, Modified NR with Armijo line search, Levenberg-Marquardt, Powell dog-leg, Broyden, Broyden II), created through the `makeBVPSolver()` factory from the `mathOp::rootFindingMethod` enum. `CTR` owns one via `unique_ptr` and recreates it on copy (solvers are stateless).
- **`boostBlazeAlgebra.hpp`** — custom `BlazeBVPAlgebra` gluing Boost.Odeint steppers to Blaze vectors. Kept in namespace `ctr` deliberately (injecting into `boost::numeric::odeint` would be UB); it is passed explicitly as a stepper template argument.
- **`Observer.hpp`/`Observer.cpp`** — vestigial. The class was removed (replaced by a lambda inside `CTR::ODESolver`); `Observer.cpp` is not in the CMake source list.

### Conventions

- Units are SI throughout: meters, Pascals, radians. Angles from degrees via `mathOp::deg2Rad`.
- Joint configuration `q` is 6-DOF: `[β₁, β₂, β₃, α₁, α₂, α₃]` — linear retractions (negative, in meters) then rotations. Tubes are ordered innermost-first everywhere (`Tube` arrays, state indices, matrices).
- Linear algebra is Blaze (`blaze::StaticVector`/`StaticMatrix`/`HybridMatrix`, mostly `columnMajor`); avoid introducing other math types at API boundaries.
- Doxygen comments on all public API; documentation is published from `docs/html` (GitHub Pages). The Doxygen CMake target is currently commented out in the root `CMakeLists.txt`.
