# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

---

[//]: # (Image References)

[image1]: ./writeup_images/groundtruth.est.sensor.track.png "Ground truth to estimate comparison"
[image2]: ./writeup_images/NIS_lidar.png "NIS lidar"
[image3]: ./writeup_images/NIS_radar.png "NIS RADAR"
[image4]: ./writeup_images/P0134.diagonal.stability.png "Stability of diagonal of P_ matrix"
[image5]: ./writeup_images/P2diagonal.stability.png "Stability of diagonal of P_ matrix P2,2"
[image6]: ./writeup_images/position_accuracy.png "Position accuracy"
[image7]: ./writeup_images/velocity_accuracy.png "Velocity accuracy"

## Rubric points

### Compilation

Code compiles without errors.

Note: it does produce a warning on a missing link directory, I believe this is
an artifact of my setup.

```
(carnd-term1) Chris-MacBook-Pro:build chisolm$ make && ./UnscentedKF | tee t 
Scanning dependencies of target UnscentedKF
[ 25%] Building CXX object CMakeFiles/UnscentedKF.dir/src/ukf.cpp.o
[ 50%] Building CXX object CMakeFiles/UnscentedKF.dir/src/main.cpp.o
[ 75%] Building CXX object CMakeFiles/UnscentedKF.dir/src/tools.cpp.o
[100%] Linking CXX executable UnscentedKF
ld: warning: directory not found for option '-L/usr/local/Cellar/libuv/1.11.0/lib'
[100%] Built target UnscentedKF
```

### Accuracy

Plot overlay of ground truth, sensor measurements and estimate position.

![alt text][image1]

Accuracy is within limits for the data set.

RMSE DataSet 1:
```
X: 0.0794
Y: 0.0843
VX: 0.3935
VY: 02449
```

Although apparently not a requirement for the rubric, the RMSE values for the second data
set in the simulator are provided below.

RMSE Dataset 2:
```
X: 0.0905
Y: 0.0811
VX: 0.5765
VY: 0.2873
```

Position accuracy was excellent.

![alt text][image6]

Velocity estimated accuracy suffers significant at the beginning of the run.

![alt text][image7]

### Sensor Fusion Algorithm

My sensor fusion algorithm is taken nearly largely from the class lessons and follows 
the same intended flow.

### First Measurements

My filter initializes position values from the first sensor data.  It does not initialize
the velocity vector as that information is not available.

### Kalman Filter Predict then Update

It follows the given order.

### Radar and Lidar measurements

It uses both types of measurements.  And obeys the use_laser and use_radar_ flags for use
of each data type.

NIS measurements for both LIDAR and RADAR:

![alt text][image2]
![alt text][image3]

### Code efficiency

I moved most/all of the static computations to the class initialization.  In some cases
matrix allocation as well, even in cases where the matrix is not used elsewhere in the
project.  There still a significant number of matrix object creations that could be
moved out of the code execution flow.

### Significant problem encountered

I had a significant problem that caused the P_ matrix to have a negative eigan value 
and the llt() transform failed in that case as it's only valid with a positive definite
matrix.

It returns a matrix that is some corrupted(partially transformed?) version of the
original matrix.

The source of that problem was that sigma points were translation into measurement space(Zsig_) the
angle was near -pi/pi with some sigma points on both sides. So a simple sum for z_pred gives a nearly random vector.  Which in turn causes an unusual S_ matrix.

```
Zsig_ i  0 5.49461 -3.1405  2.1833
Zsig_ i  1  5.28464 -3.13069  2.04544
Zsig_ i  2 5.50581 3.10736 2.10686
Zsig_ i  3  5.50156 -3.13655  2.27089
Zsig_ i  4  5.46266 -3.13813  1.84572
Zsig_ i  5 5.49205 -3.1403 2.13301
Zsig_ i  6  5.49686 -3.13957  2.23369
Zsig_ i  7 5.49461 -3.1405 2.16697
Zsig_ i  8 5.70524 3.13348 2.32565
Zsig_ i  9   5.4889 -3.10507  2.26021
Zsig_ i  10 5.48756 3.13874  2.0953
Zsig_ i  11 5.52569 3.13997 2.50907
Zsig_ i  12  5.49717 -3.14071  2.23334
Zsig_ i  13  5.49237 -3.14143  2.13309
Zsig_ i  14 5.49461 -3.1405 2.19961

z_pred_ 5.49513 -0.627191 2.18257

z_diff i  0 -0.000515722     -2.51331  0.000735618
z_diff i  1 -0.210491   -2.5035 -0.137124
z_diff i  2  0.0106825   -2.54864 -0.0757059
z_diff i  3 0.00642894   -2.50936  0.0883285
z_diff i  4 -0.0324639   -2.51094  -0.336848
z_diff i  5 -0.00308155     -2.5131  -0.0495543
z_diff i  6 0.00173391   -2.51238  0.0511283
z_diff i  7 -0.000515722     -2.51331   -0.0155982
z_diff i  8 0.210109 -2.52252 0.143082
z_diff i  9 -0.00622547    -2.47788   0.0776473
z_diff i  10 -0.00756569    -2.51726   -0.087265
z_diff i  11 0.0305601  -2.51603  0.326501
z_diff i  12 0.00204233   -2.51352  0.0507785
z_diff i  13 -0.00276058    -2.51424  -0.0494727
z_diff i  14 -0.000515722     -2.51331    0.0170449

final P_ =
  0.00412003   0.00182563   -0.0071874 -0.000761704  -0.00206853
  0.00182563    0.0032634  -0.00163036  0.000566004   0.00101482
  -0.0071874  -0.00163036    0.0122538  -0.00153668  -0.00141272
-0.000761704  0.000566004  -0.00153668 -3.01643e-05  0.000333001
 -0.00206853   0.00101482  -0.00141272  0.000333001   0.00369449
```

I ended up implementing a solution where all the sigma points are pushed to
one side of the -pi/pi boundary, even if they are slightly beyond pi.  The
common 'direction' of pi produces a useful prediction vector and does
not cause additional problems.

### Final  

My implementation works and meets the accuracy goals.

Possible improvements could be made in the early cycles.  Plotted below are 
the diagonal elements of the P_ matrix.  P(2,2) is plotted separately for scale.

![alt text][image4]
![alt text][image5]

The position and velocity components settle very early, however the yaw and yaw
rate components take 40-50 measurements to settle to the range the stay for the
remainder of the run.  Possible changes may be initialization closer to their
stable range or some method to cause a faster change in the early cycles.

Velocity estimates could still be improved in the first 0-7 measurements.



