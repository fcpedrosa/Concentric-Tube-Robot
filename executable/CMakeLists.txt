add_executable(ctr_general_use Main.cpp)
set_target_properties(ctr_general_use PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

target_link_libraries(ctr_general_use
    PRIVATE
    CTRlib
    ctr_project_options
    ctr_sanitizers
    Threads::Threads
)