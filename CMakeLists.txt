cmake_minimum_required(VERSION 3.22)
project(CTR_General_Use VERSION 2.1.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# addind the library and main executable directories to the project tree
add_subdirectory(ctr_library) # static CTR library
add_subdirectory(executable)  # main executable

# Boost
find_package(Boost REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost found")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost not found")
endif()

# BLAS
set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC ON)
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message(STATUS "BLAS found")
else()
    message(FATAL_ERROR "BLAS not found")
endif()

# LAPACK
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    message(STATUS "LAPACK found")
else()
    message(FATAL_ERROR "LAPACK not found")
endif()

# Blaze
find_package(Blaze REQUIRED)
if(Blaze_FOUND)
    message(STATUS "Blaze found")
    target_link_libraries(${PROJECT_NAME} INTERFACE blaze::blaze)
else()
    message(FATAL_ERROR "Blaze not found")
endif()

# TBB
find_package(TBB REQUIRED)
if(TBB_FOUND)
    message(STATUS "TBB found")
    target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb)
else()
    message(FATAL_ERROR "TBB not found")
endif()

# Link all necessary libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        ${Boost_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
)

# Doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

# Optionally enable CPack if needed
# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)