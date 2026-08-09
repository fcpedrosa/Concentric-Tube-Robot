add_executable(ctr_benchmarks bench_ctr.cpp)

target_link_libraries(ctr_benchmarks
    PRIVATE
    ctr::kinematics
    ctr_project_options
)

# Reuse the shared test fixture (reference robot definition).
target_include_directories(ctr_benchmarks PRIVATE ${PROJECT_SOURCE_DIR}/tests/support)
