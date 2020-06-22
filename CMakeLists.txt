cmake_minimum_required(VERSION 3.11)
project (tensor_network_qcircuit CXX)

enable_testing()

set(ITENSOR_DIR "external/itensor")
set(CMAKE_CXX_FLAGS "-std=c++17 -O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++17 -DDEBUG")
# DDEBUG is an option for ITensor.
# If not defining the flag, EVERY assert() will be disabled.

include_directories("${ITENSOR_DIR}")
include_directories("./include")
link_directories("${ITENSOR_DIR}/lib")

## Code runner configuration
add_executable(main mps_mylib_test.cpp)
target_link_libraries(main -litensor -llapack -lblas -lpthread)


## C++ test configuration
add_subdirectory(test)


## Pybind configuration
add_subdirectory("external/pybind11")
file(GLOB pybind_src "src/*.cpp")
pybind11_add_module(qcircuit ${pybind_src})
target_link_libraries(qcircuit PRIVATE -litensor -llapack -lblas -lpthread)
