# listing the CTR library source files into a variable
set(CTR_LIB_SOURCES
    src/CTR.cpp
    src/Observer.cpp
    src/ODESystem.cpp
    src/Segment.cpp
    src/Tube.cpp
)

# creating a CTR static library
add_library(CTR STATIC ${CTR_LIB_SOURCES})
target_include_directories(CTR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)