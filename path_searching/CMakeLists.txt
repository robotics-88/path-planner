cmake_minimum_required(VERSION 3.5)
project(path_searching)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

find_package(ament_cmake REQUIRED)

find_package(rclcpp)
find_package(std_msgs)
find_package(visualization_msgs)
find_package(cv_bridge)
find_package(sensor_msgs)
find_package(pcl_ros)
find_package(rcl_interfaces)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

include_directories(  
  include 
  ${EIGEN3_INCLUDE_DIRS} 
  ${PCL_INCLUDE_DIRS}
)

add_library(path_searching 
    src/kinodynamic_astar.cpp 
)
target_link_libraries(path_searching
    ${PCL_LIBRARIES}
)

ament_target_dependencies(path_searching
  rclcpp
  sensor_msgs
  pcl_conversions
  pcl_ros
  rcl_interfaces
)

ament_export_include_directories(include)
ament_export_libraries(path_searching)

install(DIRECTORY include/
  DESTINATION include
)

install(TARGETS
  path_searching
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib/${PROJECT_NAME}
)

ament_package()