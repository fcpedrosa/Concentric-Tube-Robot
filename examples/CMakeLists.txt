add_executable(ctr_demo fk_ik_demo.cpp)

target_link_libraries(ctr_demo
    PRIVATE
    ctr::kinematics
    ctr_project_options
    ctr_sanitizers
)
