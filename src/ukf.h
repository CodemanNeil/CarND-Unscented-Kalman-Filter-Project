#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
private:

    /**
     * Sets Xsig_aug_ to new augmented sigma points.  Called at beginning of Prediction().
     */
    void createAugmentedSigmaPoints();

    /**
     * Uses Xsig_aug_ to predict sigma points and set in Xsig_pred_.
     * Called following createAugmentedSigmaPoints().
     * @param delta_t time between last and current measurements in seconds
     */
    void predictSigmaPoints(long double delta_t);

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

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    ///* Augmented state vector (length: n_aug_)
    VectorXd x_aug_;

    ///* Augmented covariance matrix (n_aug_ * n_aug_)
    MatrixXd P_aug_;

    ///* Augmented sigma points matrix (n_aug_ * n_aug_)
    MatrixXd Xsig_aug_;

    ///* R measurement noise covariance matrix for radar
    MatrixXd R_radar_;

    ///* R measurement noise covariance matrix for lidar
    MatrixXd R_laser_;

    ///* H measurement matrix for LIDAR
    MatrixXd H_;

    ///* previous timestamp of last measurement, in us
    long long previous_timestamp_;

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

    ///* Number of sigma points created
    int n_sigma_;

    ///* Radar measurement dimension
    int n_z_radar_;

    ///* Laser measurement dimension
    int n_z_laser_;

    ///* Sigma point spreading parameter
    double lambda_;

    ///* the current NIS for radar
    double NIS_radar_;

    ///* the current NIS for laser
    double NIS_laser_;

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
     * prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void prediction(long double delta_t);

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void updateLidar(MeasurementPackage meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void updateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
