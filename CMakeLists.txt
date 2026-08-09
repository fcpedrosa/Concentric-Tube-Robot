cmake_minimum_required(VERSION 3.26)

project(ctr_kinematics
    VERSION 3.0.0
    DESCRIPTION "Cosserat-rod kinematics for three-tube concentric tube robots"
    HOMEPAGE_URL "https://github.com/fcpedrosa/Concentric-Tube-Robot"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ─── Options ──────────────────────────────────────────────────────────────────
option(CTR_ENABLE_WARNINGS "Enable project-wide warnings" ${PROJECT_IS_TOP_LEVEL})
option(CTR_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
option(CTR_ENABLE_SANITIZERS "Enable Address/Undefined sanitizers" OFF)
option(CTR_ENABLE_NATIVE "Build with -march=native" OFF)
option(CTR_ENABLE_LTO "Enable interprocedural optimization (LTO)" OFF)
option(CTR_BUILD_EXAMPLES "Build example executables" ${PROJECT_IS_TOP_LEVEL})
option(CTR_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(CTR_BUILD_DOCS "Build Doxygen documentation target" OFF)

if(PROJECT_IS_TOP_LEVEL)
    include(CTest) # defines BUILD_TESTING (default ON)
endif()

# Default to a Release build when nothing was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Output directories (build-tree convenience only)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ─── Dependencies ─────────────────────────────────────────────────────────────
# The library needs exactly two third-party packages:
#   Boost (header-only use: Boost.Odeint) and Blaze (linear algebra).
find_package(Boost 1.74 REQUIRED CONFIG)
find_package(blaze REQUIRED CONFIG)

# TEMPORARY: mathOp::pInv computes an SVD via blaze::svd, which calls LAPACK
# (dgesdd) unconditionally. pInv is scheduled for removal (replaced by
# fixed-size direct solves); drop this dependency together with it.
find_package(LAPACK REQUIRED)

if(CTR_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ctr_ipo_ok OUTPUT ctr_ipo_msg)
    if(ctr_ipo_ok)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "IPO/LTO not supported: ${ctr_ipo_msg}")
    endif()
endif()

# ─── Internal interface targets (never installed/exported) ────────────────────
add_library(ctr_project_options INTERFACE)

if(CTR_ENABLE_WARNINGS)
    if(MSVC)
        target_compile_options(ctr_project_options INTERFACE /W4 $<$<BOOL:${CTR_WARNINGS_AS_ERRORS}>:/WX>)
    else()
        target_compile_options(
            ctr_project_options
            INTERFACE
            -Wall
            -Wextra
            -Wpedantic
            $<$<BOOL:${CTR_WARNINGS_AS_ERRORS}>:-Werror>
        )
    endif()
endif()

if(CTR_ENABLE_NATIVE AND NOT MSVC)
    target_compile_options(ctr_project_options INTERFACE -march=native)
endif()

add_library(ctr_sanitizers INTERFACE)

if(CTR_ENABLE_SANITIZERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(ctr_sanitizers INTERFACE -fsanitize=address -fsanitize=undefined)
        target_link_options(ctr_sanitizers INTERFACE -fsanitize=address -fsanitize=undefined)
    elseif(MSVC)
        message(WARNING "Sanitizers are not currently configured for MSVC.")
    endif()
endif()

# ─── Project tree ─────────────────────────────────────────────────────────────
add_subdirectory(ctr_library)

if(CTR_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTING)
    add_subdirectory(tests)
endif()

if(CTR_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
