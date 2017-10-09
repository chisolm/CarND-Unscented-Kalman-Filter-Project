#ifndef UKF_H
#define UKF_H

#include <vector>
#include <string>
#include <fstream>
#include "measurement_package.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* aumented state vector: [x, nua, nupsi ] in SI units and rad
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_aug_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  int n_z_radar_;
  int n_z_lidar_;

  ///* Sigma point spreading parameter
  double lambda_;

  // previous timestamp
  double previous_timestamp_;

  double sqrt_lambda_n_aug;
  VectorXd weights;

  MatrixXd R_;
  MatrixXd R_lidar_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig_;
  MatrixXd Zsig_lidar_;

  // mean predicted measurement
  VectorXd z_pred_;
  VectorXd z_pred_lidar_;

  // measurement covariance matrix S
  MatrixXd S_;
  MatrixXd S_lidar_;

  double NIS_radar_;
  double NIS_lidar_;

  std::ofstream vis_out_file;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void AugmentedSigmaPoints();
  void SigmaPointPrediction(double delta_t);
  void PredictMeanAndCovariance();

  void PredictRadarMeasurement();
  void UpdateRadarState(MeasurementPackage meas_package);

  void PredictLidarMeasurement();
  void UpdateLidarState(MeasurementPackage meas_package);
};

#endif /* UKF_H */
