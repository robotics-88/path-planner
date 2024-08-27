# cmake_minimum_required(VERSION 3.5)
# project(fast_lio)

# SET(CMAKE_BUILD_TYPE "Debug")

# ADD_COMPILE_OPTIONS(-std=c++14 )
# ADD_COMPILE_OPTIONS(-std=c++14 )
# set( CMAKE_CXX_FLAGS "-std=c++14 -O3" )

# add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fexceptions" )
# set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread -std=c++0x -std=c++14 -fexceptions")

# option(DEPLOY "Deploy on hardware" OFF)
# if(DEPLOY)
#   add_definitions(-DDEPLOY)
# endif(DEPLOY)

# find_package(OpenMP QUIET)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")

# find_package(PythonLibs REQUIRED)
# find_path(MATPLOTLIB_CPP_INCLUDE_DIRS "matplotlibcpp.h")


# find_package(ament_cmake REQUIRED)

# find_package(rclcpp REQUIRED)
# find_package(fast_lio_msgs REQUIRED)
# find_package(nav_msgs REQUIRED)
# find_package(sensor_msgs REQUIRED)
# find_package(std_msgs REQUIRED)
# find_package(pcl_ros REQUIRED)
# find_package(tf2 REQUIRED)
# find_package(tf2_ros)
# find_package(tf2_eigen_kdl REQUIRED)

# find_package(Eigen3 REQUIRED)
# find_package(PCL REQUIRED)
# find_package(OpenCV REQUIRED)

# include_directories(
#   ${EIGEN3_INCLUDE_DIR}
#   ${PCL_INCLUDE_DIRS}
#   ${PYTHON_INCLUDE_DIRS}
#   include)

# add_executable(loam_feat_extract src/feature_extract.cpp)
# target_link_libraries(loam_feat_extract
#   ${PCL_LIBRARIES}
# )
# ament_target_dependencies(loam_feat_extract
#   rclcpp
#   sensor_msgs
#   pcl_ros
#   fast_lio_msgs
# )


# add_executable(loam_laserMapping src/laserMapping.cpp include/kd_tree/kd_tree.cpp)
# target_link_libraries(loam_laserMapping
#   ${PCL_LIBRARIES} 
#   ${OpenCV_LIBS} 
#   ${PYTHON_LIBRARIES})
# target_include_directories(loam_laserMapping PRIVATE ${PYTHON_INCLUDE_DIRS})
# # target_compile_definitions(loam_laserMapping PRIVATE ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

# # add_executable(kd_tree_test include/kd_tree/kd_tree.cpp src/kd_tree_test.cpp)
# # target_link_libraries(kd_tree_test ${PCL_LIBRARIES})
# ament_target_dependencies(loam_laserMapping
#   rclcpp
#   fast_lio_msgs
#   geometry_msgs
#   sensor_msgs
#   nav_msgs
#   tf2_ros
#   tf2_eigen_kdl
#   pcl_ros
# )

# install(TARGETS
#     loam_feat_extract
#     loam_laserMapping
#     DESTINATION lib/${PROJECT_NAME}
# )

# install(DIRECTORY launch DESTINATION share/${PROJECT_NAME})

# ament_package()
