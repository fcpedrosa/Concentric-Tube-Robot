add_executable(${PROJECT_NAME} Main.cpp)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        CTRlib
        project_options
        project_warnings
        pthread
        TBB::tbb
        blaze::blaze
        ${Boost_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
)

# Ensure all necessary include directories are available for the executable
target_include_directories(${PROJECT_NAME} PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${BLAS_INCLUDE_DIRS}
    ${LAPACK_INCLUDE_DIRS}
    ${BLAZE_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
)