add_library(ctr_kinematics STATIC
    src/CTR.cpp
    src/BVPSolver.cpp
    src/ODESystem.cpp
    src/Segment.cpp
    src/Tube.cpp
)
add_library(ctr::kinematics ALIAS ctr_kinematics)

set_target_properties(ctr_kinematics PROPERTIES
    VERSION ${PROJECT_VERSION}
    EXPORT_NAME kinematics)

target_include_directories(ctr_kinematics PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(ctr_kinematics PUBLIC cxx_std_23)

# Blaze's shared-memory parallelization only adds dispatch overhead (and an
# OpenMP dependency) for this library's small fixed-size kernels. Blaze is
# header-only, so the setting must be PUBLIC: consumers have to compile
# against the same configuration.
target_compile_definitions(ctr_kinematics PUBLIC
    BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)

target_link_libraries(ctr_kinematics
    PUBLIC
        Boost::headers
        blaze::blaze
    PRIVATE
        "$<BUILD_LOCAL_INTERFACE:ctr_project_options>"
        "$<BUILD_LOCAL_INTERFACE:ctr_sanitizers>")

# ─── Install & export ─────────────────────────────────────────────────────────
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS ctr_kinematics
    EXPORT ctr_kinematics-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ctr_kinematics-targets
    NAMESPACE ctr::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ctr_kinematics)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/ctr_kinematics-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ctr_kinematics-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ctr_kinematics)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ctr_kinematics-config-version.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ctr_kinematics-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ctr_kinematics-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ctr_kinematics)
