cmake_minimum_required(VERSION 3.21)
project(MatrixLib VERSION 1.0.0 LANGUAGES CXX)

option(LAB7_INJECT_BUG "Intentionally break submatrix replacement for the experiment" OFF)
option(MATRIXLIB_FETCH_GTEST "Download GoogleTest if no installed package is found" ON)

function(matrixlib_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /fp:precise)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

add_library(MatrixLib STATIC src/Matrix.cpp)
target_compile_features(MatrixLib PUBLIC cxx_std_20)
set_target_properties(MatrixLib PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(MatrixLib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
matrixlib_warnings(MatrixLib)

if(LAB7_INJECT_BUG)
    target_compile_definitions(MatrixLib PRIVATE LAB7_INJECT_BUG=1)
endif()

add_executable(MatrixDemo examples/main.cpp)
target_link_libraries(MatrixDemo PRIVATE MatrixLib)
set_target_properties(MatrixDemo PROPERTIES CXX_EXTENSIONS OFF)
matrixlib_warnings(MatrixDemo)

include(CTest)
if(BUILD_TESTING)
    find_package(GTest CONFIG QUIET)
    if(NOT TARGET GTest::gtest_main)
        if(NOT MATRIXLIB_FETCH_GTEST)
            message(FATAL_ERROR "GoogleTest not found. Install it with vcpkg or enable MATRIXLIB_FETCH_GTEST.")
        endif()
        include(FetchContent)
        if(MSVC)
            set(gtest_force_shared_crt ON CACHE BOOL "Match the MSVC runtime used by MatrixTests" FORCE)
        endif()
        FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG b514bdc898e2951020cbdca1304b75f5950d1f59 # v1.15.2
        )
        FetchContent_MakeAvailable(googletest)
    endif()

    add_executable(MatrixTests tests/MatrixTests.cpp)
    target_include_directories(MatrixTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    target_link_libraries(MatrixTests PRIVATE MatrixLib GTest::gtest_main)
    set_target_properties(MatrixTests PROPERTIES CXX_EXTENSIONS OFF)
    matrixlib_warnings(MatrixTests)

    include(GoogleTest)
    gtest_discover_tests(MatrixTests)
endif()
