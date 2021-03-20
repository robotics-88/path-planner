# dyn_small_obs_avoidance

# Avoiding dynamic small obstacles with onboard sensing and computating on aerial robots

This repository is used for UAV dynamic small obstacles avoidance. It is a complete system for lidar-based UAV, including FAST-LIO slam, time-accumulated KD-Tree mapping and kinodynamic A* search modules. It is able to avoid dynamic small obstacles (down to 20mm diameter bars) by running at 50Hz.

![RAL_coverfigure6](https://user-images.githubusercontent.com/23183555/109411264-d39ccf00-79db-11eb-85d3-764a6b17235d.png)

![outdoor2_movingobs](https://user-images.githubusercontent.com/23183555/109411282-f202ca80-79db-11eb-8ae0-704bda25ff12.png)

The source code will be available by April.

Related paper:
"Avoiding dynamic small obstacles with onboard sensing and computating on aerial robots", 
available on arxiv now https://arxiv.org/abs/2103.00406.

Related video:
https://youtu.be/pBHbQ_J1Qhc

 
## 1. Prerequisites
As the same as the prerequisites as FAST-LIO.
### 1.1 **Ubuntu** and **ROS**
Ubuntu >= 18.04.

ROS    >= Melodic. [ROS Installation](http://wiki.ros.org/ROS/Installation)

### 1.2. **PCL && Eigen**
PCL    >= 1.8,   Follow [PCL Installation](http://www.pointclouds.org/downloads/linux.html).

Eigen  >= 3.3.4, Follow [Eigen Installation](http://eigen.tuxfamily.org/index.php?title=Main_Page).

### 1.3. **livox_ros_driver**
Follow [livox_ros_driver Installation](https://github.com/Livox-SDK/livox_ros_driver).


## 2. Build
Clone the repository and catkin_make:

```
    cd ~/catkin_ws/src
    git clone https://github.com/hku-mars/dyn_small_obs_avoidance.git
    cd ..
    catkin_make
    source devel/setup.bash
```
# 3.Run demo
### 3.１ Start program
```
    source devel/setup.bash
    roslaunch path_planning demo.launch 
```

### 3.2 Run rosbag or directly fly
You can download our demos rosbag from https://drive.google.com/drive/folders/1knQwnrbwunGIvXzOL6wWKkCtOkGE22bG?usp=sharing
And play the bag by:
```
    rosbag play XXX.bag
```

### 3.3 Set target point
```
    rostopic pub /goal geometry_msgs/PoseStamped '{header: {stamp: now, frame_id: "camera_init"}, pose: {position: {x: 5.0, y: 0.0, z: 1.0}, orientation: {w: 1.0}}}'

```

## 5.Acknowledgments
Thanks for FAST-PLANNER(Zhou, Boyu and Gao, Fei and Wang, Luqi and Liu, Chuhao and Shen, Shaojie. Robust and efficient quadrotor trajectory generation for fast autonomous flight), [FAST-PLANNER](https://github.com/HKUST-Aerial-Robotics/Fast-Planner.git).