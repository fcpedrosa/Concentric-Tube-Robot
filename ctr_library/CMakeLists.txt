# Define project name and version
project(CTR VERSION 1.0 LANGUAGES CXX)

# listing the CTR library source files into a variable
set(CTR_LIB_SOURCES
    src/CTR.cpp
    src/Observer.cpp
    src/ODESystem.cpp
    src/Segment.cpp
    src/Tube.cpp
)

# Create the CTR static library
add_library(CTRlib STATIC ${CTR_LIB_SOURCES})

# Include directories for the CTR library
target_include_directories(CTRlib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)