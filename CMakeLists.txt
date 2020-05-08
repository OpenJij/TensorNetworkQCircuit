cmake_minimum_required(VERSION 3.11)
project (tensor_network_qcircuit CXX)

set(ITENSOR_DIR "external/itensor")
set(CMAKE_CXX_FLAGS "-std=c++17 -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++17 -DDEBUG")

include_directories("${ITENSOR_DIR}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g -O0")
link_directories("${ITENSOR_DIR}/lib")

add_executable(main mps_mylib_test.cc)
target_link_libraries(main -lblas -llapack -lpthread -litensor)
