#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	// TODO need to tune
	std_a_ = 0.5;  // 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	// TODO need to tune
	std_yawdd_ = 0.3;  // 30;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	previous_timestamp_ = 0;

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	///* Sigma point spreading parameter
	lambda_ = 3 - n_x_;

	//create example sigma point matrix
	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//create augmented mean vector
	x_aug_ = VectorXd(n_aug_);

	//create augmented state covariance
	P_aug_ = MatrixXd(n_aug_, n_aug_);

	// Precomputed values...
	sqrt_lambda_n_aug = sqrt(lambda_+n_aug_);

	n_z_radar_ = 3;
	n_z_lidar_ = 2;
	//create matrix for sigma points in measurement space
	Zsig_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
	Zsig_lidar_ = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

	//mean predicted measurement
	z_pred_ = VectorXd(n_z_radar_);
	z_pred_lidar_ = VectorXd(n_z_lidar_);

	//measurement covariance matrix S
	S_ = MatrixXd(n_z_radar_,n_z_radar_);
	S_lidar_ = MatrixXd(n_z_lidar_,n_z_lidar_);


	R_ = MatrixXd(n_z_radar_,n_z_radar_);
	R_ << std_radr_*std_radr_, 0, 0,
				0, std_radphi_*std_radphi_, 0,
				0, 0, std_radrd_*std_radrd_;
	R_lidar_ = MatrixXd(n_z_lidar_,n_z_lidar_);
	R_lidar_ << std_laspx_*std_laspx_, 0,
				0, std_laspy_*std_laspy_;

	//create vector for weights
	weights = VectorXd(2*n_aug_+1);
	weights.fill(0.5/(n_aug_+lambda_));
	weights(0) = lambda_/(lambda_+n_aug_);

	NIS_lidar_ = 0.0;
	NIS_radar_ = 0.0;


	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	/*
	*/
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		cout << "RADAR " << meas_package.raw_measurements_[0] << ", " <<
							meas_package.raw_measurements_[1] << ", " <<
							meas_package.raw_measurements_[2] << endl;
	} else {
		// Lidar updates
		cout << "LASER " << meas_package.raw_measurements_[0] << ", " <<
							meas_package.raw_measurements_[1] << endl;
	}

	if (!is_initialized_) {

		P_ << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0,
			0, 0, .1, 0, 0,
			0, 0, 0, .1, 0,
			0, 0, 0, 0, .1;

        MatrixXd Pxx;
        Pxx = MatrixXd(7, 7);
        Pxx << -0.00372461, 0.0539572, -0.302925, -0.0701267, -0.471767, 0, 0,
                0.0539572, -0.100754, 0.240057, 0.171488, 0.997984, 0, 0,
                -0.302925, 0.240057, 5.66721, -0.718579, -3.62293, 0, 0,
                -0.0701267, 0.171488, -0.718579, -0.212282, -1.05473, 0, 0,
                -0.471767, 0.997984, -3.62293, -1.05473, -1.02394, 0, 0,
                0, 0, 0, 0, 0, 900, 0,
                0, 0, 0, 0, 0, 0, 900;

        std::cout << "Pxx " << std::endl << Pxx << std::endl;
        //create square root matrix
        MatrixXd A = Pxx.llt().matrixL();
        std::cout << "A " << std::endl << A << std::endl;
		cout << "To check this, let us compute A * A.transpose()" << endl;
		cout << A * A.transpose() << endl;

		Pxx << 0.0188397, -0.00087777, -0.00321798, -0.00944216, -0.0503088, 0, 0,
			-0.00087777, 0.0151152, 0.0110446, 0.0254677, 0.120575, 0, 0,
			-0.00321798, 0.0110446, 0.757013, 0.0449782, 0.539508, 0, 0,
			-0.00944216, 0.0254677, 0.0449782, 0.134377, 1.1674, 0, 0,
			-0.0503088, 0.120575, 0.539508, 1.1674, 14.4438, 0, 0,
        	0, 0, 0, 0, 0, 900, 0,
        	0, 0, 0, 0, 0, 0, 900;

        std::cout << "Pxx " << std::endl << Pxx << std::endl;
        //create square root matrix
        A = Pxx.llt().matrixL();
        std::cout << "A " << std::endl << A << std::endl;
		cout << "To check this, let us compute A * A.transpose()" << endl;
		cout << A * A.transpose() << endl;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    	/**
    	Convert radar from polar to cartesian coordinates and initialize state.
    	*/
    	x_ << meas_package.raw_measurements_[0] *
                  cos(meas_package.raw_measurements_[1]),
                  meas_package.raw_measurements_[0] *
                  sin(meas_package.raw_measurements_[1]),
                  0, 0, 0;
    	} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    	/**
    	Initialize state.
    	*/
    	x_ << meas_package.raw_measurements_[0],
                meas_package.raw_measurements_[1], 0, 0, 0;
    	}
    	previous_timestamp_ = meas_package.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == false) {
		return;
	}
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == false) {
		return;
	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	// dt - expressed in seconds
	double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = meas_package.timestamp_;
	if (dt > 30 || dt < .001) {
			cout << "Time stamp " << dt << 
						"seconds, does not make sense.  exiting." << endl;
			exit(-1);
	}

	// TODO more initialization?
	cout << "Delta T = " << endl << dt << endl;
	cout << "Before Predition x_ = " << x_.transpose() << endl;
	cout << "Before Predition P_ = " << endl << P_ << endl;

	Prediction(dt);
	cout << "After Predition x_ = " << x_.transpose() << endl;
	cout << "After Predition P_ = " << endl << P_ << endl;

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		// cout << "RADAR " << z[0] << ", " << z[1] << ", " << z[2] << endl;
		// cout << "Hj " << Hj_ << endl;
		cout << "Update RADAR " << meas_package.raw_measurements_[0] << ", " <<
							meas_package.raw_measurements_[1] << ", " <<
							meas_package.raw_measurements_[2] << endl;
		UpdateRadar(meas_package);

	} else {
		// Lidar updates
		cout << "Update LASER " << meas_package.raw_measurements_[0] << ", " <<
							meas_package.raw_measurements_[1] << endl;
		UpdateLidar(meas_package);
	}

	// print the output
	cout << "final x_ = " << endl << x_ << endl;
	cout << "final P_ = " << endl << P_ << endl;


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/

	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
	PredictLidarMeasurement();
	UpdateLidarState(meas_package);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::PredictLidarMeasurement() {

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		// double v = Xsig_pred_(2, i);
		// double yaw = Xsig_pred_(3,i);
		// double yawd = Xsig_pred_(4,i);

		Zsig_lidar_(0, i) = px;
		Zsig_lidar_(1, i) = py;
	}

	//calculate mean predicted measurement
	z_pred_lidar_.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_lidar_ = z_pred_lidar_ + weights(i) * Zsig_lidar_.col(i);
	}

	std::cout << Zsig_lidar_.rows() << "x" << Zsig_lidar_.cols() << std::endl;
	std::cout << z_pred_lidar_.rows() << "x" << z_pred_lidar_.cols() << std::endl;

	//calculate measurement covariance matrix S
	S_lidar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;

		S_lidar_ = S_lidar_ + weights(i) * z_diff * z_diff.transpose() ;
	}

	S_lidar_ = S_lidar_ + R_lidar_;
}

void UKF::UpdateLidarState(MeasurementPackage meas_package) {

	VectorXd z = VectorXd(n_z_lidar_);
	z << meas_package.raw_measurements_[0],
		meas_package.raw_measurements_[1];

	 //create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {	//iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;

		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose() ;
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z_lidar_);
	K = Tc * S_lidar_.inverse();
	VectorXd z_diff = (z - z_pred_lidar_);

	//calculate NIS
	NIS_lidar_ = z_diff.transpose() * S_lidar_.inverse() * z_diff;

	//update state mean and covariance matrix

	x_ = x_ + K * (z_diff);
	P_ = P_ - K * S_lidar_ * K.transpose();

}


void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/
	  //transform sigma points into measurement space

	PredictRadarMeasurement();
	UpdateRadarState(meas_package);
}

void UKF::PredictRadarMeasurement() {

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3,i);
		double yawd = Xsig_pred_(4,i);

		double r = sqrt(px*px + py*py);
		double phi = atan2(py, px);
		double rdot = (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt(px*px + py*py);

		Zsig_(0, i) = r;
		Zsig_(1, i) = phi;
		Zsig_(2, i) = rdot;
		std::cout << "Zsig_ i  " << i << " " << Zsig_.col(i).transpose() << std::endl;
	}

    // put aall Zsig on same side of -pi/pi border
	if (fabs(Zsig_(1, 0)) > 3.0) {
		for (int i = 1; i < 2 * n_aug_ + 1; i++) {
			std::cout << "Zsig_ i tval " << i << " " << (Zsig_(1, 0) < 0.0) << " " << (Zsig_(1, i) > 0.0) << std::endl;
			std::cout << "Zsig_ i tval " << i << " " << (Zsig_(1, 0) > 0.0) << " " << (Zsig_(1, i) < 0.0) << std::endl;
			if ((Zsig_(1, 0) < 0.0) && (Zsig_(1, i) > 0.0)) {
				Zsig_(1, i)-=2.*M_PI;
			std::cout << "Zsig_ i fix " << i << " " << Zsig_.col(i).transpose() << std::endl;
			}
			if ((Zsig_(1, 0) > 0.0) && (Zsig_(1, i) < 0.0)) {
				Zsig_(1, i)+=2.*M_PI;
			std::cout << "Zsig_ i fix " << i << " " << Zsig_.col(i).transpose() << std::endl;
			}
			std::cout << "Zsig_ i sum " << i << " " << Zsig_.col(i).transpose() << std::endl;
		}
	}

	//calculate mean predicted measurement
	z_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_ = z_pred_ + weights(i) * Zsig_.col(i);
	}

	std::cout << Zsig_.rows() << "x" << Zsig_.cols() << std::endl;
	std::cout << z_pred_.rows() << "x" << z_pred_.cols() << std::endl;
	std::cout << "z_pred_ " << std::endl << z_pred_ << std::endl;

	//calculate measurement covariance matrix S
	S_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig_.col(i) - z_pred_;


		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
		std::cout << "z_diff i  " << i << " " << z_diff.transpose() << std::endl;

		S_ = S_ + weights(i) * z_diff * z_diff.transpose() ;
	}

	std::cout << "S_  " << std::endl << S_ << std::endl;

	S_ = S_ + R_;
}

void UKF::UpdateRadarState(MeasurementPackage meas_package) {

	VectorXd z = VectorXd(3);
	z << meas_package.raw_measurements_[0],
		meas_package.raw_measurements_[1],
		meas_package.raw_measurements_[2];

	 //create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {	//iterate over sigma points

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		VectorXd z_diff = Zsig_.col(i) - z_pred_;
		//angle normalization

// TODO Again here, i need to make sure Zsig - zpred are on the same side of -pi/pi
// should be ok, adjusting Zsig above
// TODO probably x_diff(3) too
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		std::cout << "z_diff 2nd i  " << i << " " << z_diff.transpose() << std::endl;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose() ;
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z_radar_);
	K = Tc * S_.inverse();

	//update state mean and covariance matrix

	std::cout << "K " << std::endl << K << std::endl;
	std::cout << "z - z_pred_ " << std::endl << z - z_pred_ << std::endl;

	VectorXd z_diff = z - z_pred_;
// TODO Again here, i need to make sure z - zpred are on the same side of -pi/pi
	std::cout << "z_diff 3rd " << z_diff.transpose() << std::endl;
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	std::cout << "z_diff 3rd " << z_diff.transpose() << std::endl;
	NIS_radar_ = z_diff.transpose() * S_.inverse() * z_diff;

	x_ = x_ + K * (z_diff);
	P_ = P_ - K * S_ * K.transpose();

}


void UKF::AugmentedSigmaPoints() {

	// x_aug << x, std_a, std_yawdd;
	std::cout << "x_" << x_ << std::endl;
	std::cout << "x_ rows" << x_.rows() << std::endl;
	std::cout << "x_aug_ rows" << x_aug_.rows() << std::endl;
	x_aug_ << x_, 0, 0;
	std::cout << "x_aug_" << x_aug_ << std::endl;
	//std::cout << "x_aug" << std::endl << x_aug << std::endl;

	//create augmented covariance matrix
	MatrixXd P_proc_noise = MatrixXd(2, 2);
	P_proc_noise << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;

	std::cout << "P_ " << P_.rows() << "x" << P_.cols() << std::endl;
	std::cout << "P_proc_noise " << P_proc_noise.rows() << "x" << P_proc_noise.cols() << std::endl;
	std::cout << "P_aug_ " << P_aug_.rows() << "x" << P_aug_.cols() << std::endl;
	std::cout << "5/2 " << (MatrixXd::Zero(5,2)).rows() << "x" << (MatrixXd::Zero(5,2)).cols() << std::endl;
	std::cout << "P_ " << std::endl << P_ << std::endl;
	P_aug_ << P_,
				MatrixXd::Zero(5,2),
				MatrixXd::Zero(2,5),
				P_proc_noise;
	std::cout << "P_aug_ " << std::endl << P_aug_ << std::endl;
	//create square root matrix
	MatrixXd A = P_aug_.llt().matrixL();
	std::cout << "A " << std::endl << A << std::endl;
	if (A(6,6) == P_aug_(6,6)) {
		std::cout << "This is the llt() error case." << std::endl;
	}
	
	//create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;
	//std::cout << "Xsig_aug.col(0) " << Xsig_aug.col(0) << std::endl;

	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug_.col(i+1)     = x_aug_ + sqrt_lambda_n_aug * A.col(i);
	    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt_lambda_n_aug * A.col(i);
	}
}

void UKF::SigmaPointPrediction(double delta_t) {

	// Input is Xsig_aug, and output is Xsig_pred.

	//predict sigma points
	//avoid division by zero
	//write predicted sigma points into right column
	for(int i=0; i<2 * n_aug_ + 1; i++) {
		std::cout << "Xsig_aug_ i  " << i << " " << Xsig_aug_.col(i).transpose() << std::endl;
        double p_x = Xsig_aug_(0,i);
        double p_y = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        } else {
            px_p = p_x + v * cos(yaw) * delta_t;
            py_p = p_y + v * sin(yaw) * delta_t;
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
	}
}

void UKF::PredictMeanAndCovariance() {

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights(i) * Xsig_pred_.col(i);
		std::cout << "Xsig_pred_ i  " << i << " " << Xsig_pred_.col(i).transpose() << std::endl;

	}
	// std::cout << "Xsig_pred_ " << std::endl << Xsig_pred_ << std::endl;

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		std::cout << "xdiff i  " << i << " " << x_diff.transpose() << std::endl;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
	}
}