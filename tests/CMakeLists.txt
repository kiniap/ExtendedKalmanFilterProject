find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (${ExtendedKF_SOURCE_DIR}/src) 

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

add_library(Tracking tracking.cpp tracking.h)
add_executable (unitTestEKF unitTestEKF.cpp) 

target_link_libraries(unitTestEKF
                      Tracking
                      fusionEKF
                      KF
                      Tools
                      ${Boost_FILESYSTEM_LIBRARY}
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                      )

add_test(unitTestEKF unitTestEKF)