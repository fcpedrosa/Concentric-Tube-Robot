add_executable(${PROJECT_NAME} Main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC CTRlib pthread tbb)

# Ensure all necessary include directories are available for the executable
target_include_directories(${PROJECT_NAME} PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${BLAS_INCLUDE_DIRS}
    ${LAPACK_INCLUDE_DIRS}
    ${pugixml_INCLUDE_DIRS}
    ${ASSIMP_INCLUDE_DIRS}
    ${BLAZE_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
)