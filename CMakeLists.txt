cmake_minimum_required(VERSION 3.26)
project(CTR_General_Use VERSION 2.1.4 LANGUAGES CXX)

# Set the CMP0167 policy to find the upstream BoostConfig.cmake directly
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(CTR_ENABLE_WARNINGS "Enable project-wide warnings" ON)
option(CTR_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
option(CTR_ENABLE_SANITIZERS "Enable Address/Undefined sanitizers" OFF)

include(CTest)

# Set the default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set RPATH to avoid library conflicts
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Dependency discovery
find_package(Threads REQUIRED)
find_package(Boost REQUIRED)
find_package(Blaze REQUIRED)

set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC ON)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(TBB REQUIRED)

set(CTR_THIRD_PARTY_LIBS
    Boost::boost
    blaze::blaze
    TBB::tbb
    BLAS::BLAS
    LAPACK::LAPACK
    CACHE INTERNAL "List of third-party targets used across the project")

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

add_library(ctr_sanitizers INTERFACE)

if(CTR_ENABLE_SANITIZERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(ctr_sanitizers INTERFACE -fsanitize=address -fsanitize=undefined)
        target_link_options(ctr_sanitizers INTERFACE -fsanitize=address -fsanitize=undefined)
    elseif(MSVC)
        message(WARNING "Sanitizers are not currently configured for MSVC.")
    endif()
endif()

# addind the library and main executable directories to the project tree
add_subdirectory(ctr_library) # static CTR library
add_subdirectory(executable) # main executable

# # Doxygen documentation
# find_package(Doxygen)
# if(DOXYGEN_FOUND)
# set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
# set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

# add_custom_target(doc_doxygen ALL
# COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
# WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
# COMMENT "Generating API documentation with Doxygen"
# VERBATIM)
# endif()

# Optionally enable CPack if needed
# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)