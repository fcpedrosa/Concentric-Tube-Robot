cmake_minimum_required(VERSION 3.22)
project(CTR_General_Use VERSION 2.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Generated binary files will be directed to bin folder 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# addind the library and main executable directories to the project tree
add_subdirectory(ctr_library) # static CTR library
add_subdirectory(executable)  # main executable

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost REQUIRED)
find_package(LAPACK REQUIRED)

set(BLA_VENDOR OpenBLAS)
set(BLA_STATIC                   ON)  # use static linkage
find_package(BLAS)
find_package(Blaze REQUIRED)

# Check for the Boost libray, if found print message, include dirs and link libraries.
if (Boost_FOUND)
    message(STATUS "Boost has been found!")
    include_directories(${Boost_INCLUDE_DIRS})
    link_libraries(${Boost_LIBRARIES})
elseif(NOT Boost_FOUND)
    error(STATUS "Boost has not been found!")
endif()

if (BLAS_FOUND)
    message(STATUS "BLAS has been found!")
    include_directories(${BLAS_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${BLAS_LIBRARIES})
elseif(NOT BLAS_FOUND)
    error(STATUS "BLAS library has not been found")
endif()

# Check for the LAPACK libray, if found print message, include dirs and link libraries.
if (LAPACK_FOUND)
    message("LAPACK has been found!")
    include_directories(${LAPACK_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
elseif(NOT LAPACK_FOUND)
    error("LAPACK API has not been found")
endif()

# Check for the Blaze libray, if found print message, include dirs and link libraries.
find_package( blaze )
if( blaze_FOUND )
   message(STATUS "Blaze has been found!")
   target_link_libraries(${PROJECT_NAME} INTERFACE blaze::blaze)
elseif(NOT Blaze_FOUND)
    error(STATUS "Blaze has not been found!")
endif()

# Check for the TBB libray. If found, print message, include dirs, and link libraries.
find_package(TBB REQUIRED)
if( TBB_FOUND )
   message(STATUS "TBB has been found!")
   target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb)
elseif(NOT TBB_FOUND)
    error(STATUS "TBB has not been found!")
endif()

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)