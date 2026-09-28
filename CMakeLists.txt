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

# Build type and output directories are the parent project's business when
# this project is vendored via add_subdirectory().
if(PROJECT_IS_TOP_LEVEL)
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
endif()

# Must precede add_subdirectory(ctr_library): the target's INSTALL_INTERFACE
# include directory is taken from CMAKE_INSTALL_INCLUDEDIR.
include(GNUInstallDirs)

# ─── Dependencies ─────────────────────────────────────────────────────────────
# The library needs exactly two third-party packages:
#   Boost (header-only use: Boost.Odeint) and Blaze (linear algebra).
find_package(Boost 1.74 REQUIRED CONFIG)

# Blaze: prefer an installed package; otherwise fetch the 3.8.2 release (many
# distributions, e.g. Debian/Ubuntu, do not package Blaze). A fetched Blaze is
# installed alongside this library so the exported package stays consumable.
include(FetchContent)
FetchContent_Declare(blaze
    URL https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.8.2.tar.gz
    URL_HASH SHA256=4c4e1915971efbedab95790e4c5cf017d8448057fa8f8c62c46e1643bf72cbb1
    SYSTEM
    FIND_PACKAGE_ARGS 3.8 CONFIG)
# Only consulted if Blaze is fetched. Blaze's cmake_minimum_required predates
# CMP0077, so its options must be pre-seeded as cache entries: no LAPACK and
# no OpenMP (this library needs neither).
set(USE_LAPACK OFF CACHE BOOL "Blaze: require LAPACK")
set(BLAZE_SHARED_MEMORY_PARALLELIZATION OFF CACHE BOOL "Blaze: shared-memory parallelization")
set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON) # Blaze calls export(PACKAGE)
FetchContent_MakeAvailable(blaze)
if(NOT TARGET blaze::blaze)
    add_library(blaze::blaze ALIAS blaze) # Blaze's own CMake defines no alias
endif()

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

if(CTR_BUILD_DOCS)
    add_subdirectory(docs)
endif()
