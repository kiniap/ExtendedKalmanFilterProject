project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
set(CMAKE_BUILD_TYPE Debug)

set(sources src/main.cpp src/tools.cpp src/FusionEKF.cpp src/kalman_filter.cpp src/tools.h src/FusionEKF.h src/kalman_filter.h)

enable_testing()
add_subdirectory (tests)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


add_executable(ExtendedKF src/main.cpp)
add_library(fusionEKF src/FusionEKF.cpp src/FusionEKF.h)
add_library(Tools src/tools.cpp src/tools.h)
add_library(KF src/kalman_filter.cpp src/kalman_filter.h)

target_link_libraries(ExtendedKF fusionEKF Tools KF z ssl uv uWS)
