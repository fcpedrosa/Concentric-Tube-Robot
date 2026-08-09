include(FetchContent)

# Catch2 v3: prefer a system installation, fall back to FetchContent.
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
    FIND_PACKAGE_ARGS 3)
FetchContent_MakeAvailable(Catch2)

# catch_discover_tests lives in Catch.cmake (shipped in extras/ of the source
# tree; on the install side it sits next to the package config, which Catch2's
# config script adds to CMAKE_MODULE_PATH).
if(catch2_SOURCE_DIR)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
include(Catch)

# ── Unit tests: fast, no ODE integration beyond single FK solves ─────────────
add_executable(ctr_unit_tests
    unit/test_tube.cpp
    unit/test_segment.cpp
    unit/test_math.cpp
    unit/test_ctr_semantics.cpp
)
target_link_libraries(ctr_unit_tests PRIVATE ctr::kinematics Catch2::Catch2WithMain ctr_project_options ctr_sanitizers)
target_include_directories(ctr_unit_tests PRIVATE support)
catch_discover_tests(ctr_unit_tests PROPERTIES LABELS "unit")

# ── Regression tests: analytic anchors, golden values, IK, solvers ───────────
add_executable(ctr_regression_tests
    regression/test_fk_anchors.cpp
    regression/test_fk_golden.cpp
    regression/test_ik.cpp
    regression/test_solvers.cpp
)
target_link_libraries(ctr_regression_tests PRIVATE ctr::kinematics Catch2::Catch2WithMain ctr_project_options
                                                   ctr_sanitizers)
target_include_directories(ctr_regression_tests PRIVATE support)
catch_discover_tests(ctr_regression_tests PROPERTIES LABELS "regression" TIMEOUT 300)

# ── Slow tests: convergence studies etc. (excluded from quick runs via -LE slow)
add_executable(ctr_slow_tests
    slow/test_convergence.cpp
)
target_link_libraries(ctr_slow_tests PRIVATE ctr::kinematics Catch2::Catch2WithMain ctr_project_options
                                             ctr_sanitizers)
target_include_directories(ctr_slow_tests PRIVATE support)
catch_discover_tests(ctr_slow_tests PROPERTIES LABELS "slow" TIMEOUT 900)
