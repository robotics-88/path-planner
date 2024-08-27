cmake_minimum_required(VERSION 3.5)
project(path_planning)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

find_package(Eigen3 REQUIRED)
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(pcl_conversions REQUIRED)
find_package(pcl_ros REQUIRED)
find_package(std_msgs REQUIRED)
find_package(path_searching REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(trajectory_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)

find_package(pcl_ros REQUIRED) #TODO fix hack

include_directories(
  ${PROJECT_SOURCE_DIR}/../path_searching/include # TODO fix this hack

  include
  ${PCL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS} 
)

add_library(${PROJECT_NAME} SHARED
  src/path_planning.cpp
)

ament_target_dependencies(${PROJECT_NAME}
  rclcpp
  visualization_msgs
  trajectory_msgs
  nav_msgs
  sensor_msgs
  geometry_msgs
  path_searching
  pcl_ros # part of hack
)

add_executable(${PROJECT_NAME}_node src/path_planning.cpp)

target_link_libraries(${PROJECT_NAME}_node ${PROJECT_NAME})

install(TARGETS
  ${PROJECT_NAME}_node
  DESTINATION lib/${PROJECT_NAME}
)

ament_package()
