#include "ukf.h"
#include "tools.h"
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
    // Initially false
    is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // Number of states
    n_x_ = 5;

    // Number of states including augmented states
    n_aug_ = 7;

    n_sigma_ = 2*n_aug_ + 1;

    // Lambda for sigma point calculation
    lambda_ = 3 - n_aug_;

    // Radar measurement dimension
    n_z_radar_ = 3;

    // Laser measurement dimension
    n_z_laser_ = 2;

    // Initial state vector
    x_ = VectorXd(n_x_);

    // Initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    // Predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

    // Augmented state vector
    x_aug_ = VectorXd(n_aug_);

    // Augmented covariance matrix
    P_aug_ = MatrixXd(n_aug_, n_aug_);

    // Augmented sigma points matrix
    Xsig_aug_ = MatrixXd(n_aug_, n_sigma_);

    // R measurement noise covariance matrix for radar
    R_radar_ = MatrixXd::Zero(n_z_radar_, n_z_radar_);

    // R measurement noise covariance matrix for lidar
    R_laser_ = MatrixXd::Zero(n_z_laser_, n_z_laser_);

    // H measurement matrix for LIDAR
    H_ = MatrixXd(n_z_laser_, n_x_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.6;

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

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */
    weights_ = VectorXd(n_sigma_);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // Initialize variables
    if (!is_initialized_) {

        // Set initial P covariance matrix
        P_ = MatrixXd::Identity(n_x_, n_x_);

        // Set initial x state vector
        x_ << 0, 0, 0, 0, 0;

        H_ <<   1, 0, 0, 0, 0,
                0, 1, 0, 0, 0;

        // Set weights
        weights_[0] = lambda_/(lambda_ + n_aug_);
        weights_.tail(2 * n_aug_).setConstant(1./(2.*(lambda_ + n_aug_)));

        // Set R noise matrices
        R_radar_(0,0) = pow(std_radr_, 2);
        R_radar_(1,1) = pow(std_radphi_, 2);
        R_radar_(2,2) = pow(std_radrd_, 2);

        R_laser_(0,0) = pow(std_laspx_, 2);
        R_laser_(1,1) = pow(std_laspy_, 2);

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double ro = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            double ro_dot = meas_package.raw_measurements_[2];
            x_[0] = ro * cos(phi);
            x_[1] = ro * sin(phi);
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_[0] = meas_package.raw_measurements_[0];
            x_[1] = meas_package.raw_measurements_[1];
        }

        // If x,y == 0, set to very small values
        x_[0] = (x_[0] == 0) ? 1e-4 : x_[0];
        x_[1] = (x_[1] == 0) ? 1e-4 : x_[1];

        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
    } else {

        //compute the time elapsed between the current and previous measurements
        long double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;    // Convert micro_sec -> sec
        previous_timestamp_ = meas_package.timestamp_;

        /** TODO: Predict **/
        prediction(delta_t);

        /** TODO: Update **/
        if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER && use_laser_) {
            updateLidar(meas_package);
        } else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && use_radar_) {
            updateRadar(meas_package);
        }
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {long double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::prediction(long double delta_t) {
    /**
     Estimate the object's location. Modify the state vector, x_.
     Predict sigma points, the state, and the state covariance matrix.
    */

    // Generate Augmented Sigma Points
    createAugmentedSigmaPoints();

    // Predict Sigma Points
    predictSigmaPoints(delta_t);

    // Predict state mean
    x_.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
        x_ += weights_[i] * Xsig_pred_.col(i);
    }

    // Predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff[3] = Tools::normalizeAngle(x_diff[3]);
        P_ += weights_[i] * x_diff * x_diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::updateLidar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */

    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_diff = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_laser_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
    x_ = x_ + K * z_diff;
    P_ = (I - K * H_) * P_;

    // Set current laser NIS
    NIS_laser_ = Tools::calculateNIS(S, z_diff);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::updateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */

    VectorXd z = meas_package.raw_measurements_;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_radar_, n_sigma_);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_radar_);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z_radar_, n_z_radar_);

    //transform sigma points into measurement space
    for (int i = 0; i < n_sigma_; i++) {
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double psi = Xsig_pred_(3,i);

        Zsig(0,i) = sqrt(px*px + py*py);
        if (sqrt(px*px + py*py) == 0) {
            Zsig(1,i) = 0;
        } else {
            Zsig(1, i) = atan2(py, px);
        }
        Zsig(2,i) = (px*cos(psi)*v + py*sin(psi)*v)/(sqrt(px*px + py*py));
    }

    //calculate mean predicted measurement
    z_pred = Zsig * weights_;

    //calculate measurement covariance matrix S
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += (weights_[i]*z_diff*z_diff.transpose());
    }
    S += R_radar_;

    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_radar_);

    //calculate cross correlation matrix
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd x_d = Xsig_pred_.col(i) - x_;
        VectorXd z_d = Zsig.col(i) - z_pred;
        Tc += weights_[i]*(x_d*z_d.transpose());
    }
    //calculate Kalman gain K;
    MatrixXd K = Tc*S.inverse();
    VectorXd z_diff = z - z_pred;

    //update state mean and covariance matrix
    x_ += K*z_diff;
    P_ -= K*S*K.transpose();

    // Set current radar NIS
    NIS_radar_ = Tools::calculateNIS(S, z_diff);
}

/**
 * Sets Xsig_aug_ to new augmented sigma points.  Called at beginning of Prediction()
 */
void UKF::createAugmentedSigmaPoints() {
    //create augmented mean vector
    x_aug_.head(n_x_) = x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;

    //create augmented state covariance
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_,n_x_) = P_;
    P_aug_(5,5) = std_a_*std_a_;
    P_aug_(6,6) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd A = P_aug_.llt().matrixL();

    //create augmented sigma points
    Xsig_aug_.col(0)  = x_aug_;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * A.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * A.col(i);
    }
}


/**
 * Uses Xsig_aug_ to predict sigma points and set in Xsig_pred_.
 * Called following createAugmentedSigmaPoints().
 * @param delta_t time between last and current measurements in seconds
 */
void UKF::predictSigmaPoints(long double delta_t) {
    for (int i = 0; i < n_sigma_; i++) {
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yaw_dot = Xsig_aug_(4,i);
        double noise_a = Xsig_aug_(5,i);
        double noise_yaw = Xsig_aug_(6,i);
        const long long delta_t_2 = pow(delta_t, 2);

        // Create process noise vector to add to sigma point
        VectorXd xsig_noise = VectorXd(n_x_);
        xsig_noise(0) = (1./2.)*(delta_t_2)*cos(yaw)*noise_a;
        xsig_noise(1) = (1./2.)*(delta_t_2)*sin(yaw)*noise_a;
        xsig_noise(2) = delta_t*noise_a;
        xsig_noise(3) = (1./2.)*(delta_t_2)*noise_yaw;
        xsig_noise(4) = delta_t*noise_yaw;

        // Create process vector to add to sigma point
        VectorXd xsig_process = VectorXd::Zero(n_x_);
        // If divide by zero
        if (fabs(yaw_dot) > 0.001) {
            xsig_process(0) = (v/yaw_dot)*(sin(yaw + yaw_dot*delta_t) - sin(yaw));
            xsig_process(1) = (v/yaw_dot)*(-cos(yaw + yaw_dot*delta_t) + cos(yaw));
            xsig_process(3) = yaw_dot*delta_t;
        } else {
            xsig_process(0) = v*cos(yaw)*delta_t;
            xsig_process(1) = v*sin(yaw)*delta_t;
            xsig_process(3) = yaw_dot*delta_t;
        }
        // Set prediction sigma point to augmented sigma point[0-4] + process motion + proces noise
        Xsig_pred_.col(i) = Xsig_aug_.col(i).head(n_x_) + xsig_process + xsig_noise;
    }
}