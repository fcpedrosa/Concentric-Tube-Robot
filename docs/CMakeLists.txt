find_package(Doxygen REQUIRED)

add_custom_target(ctr_docs
    COMMAND Doxygen::doxygen ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen (docs/generated/html)"
    VERBATIM
)
