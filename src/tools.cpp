#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::calculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse = VectorXd::Zero(4);

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        throw std::runtime_error("Ground Truth size doesn't match estimation size");
    }

    //accumulate squared residuals
    for (unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

/**
     * Calculates NIS
     * @param S measurement covariance matrix
     * @param z_diff difference between predicted measurement and measurement vector
     * @return NIS value
     */
double Tools::calculateNIS(const MatrixXd& S, const VectorXd& z_diff) {
    return z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Normalizes angle in radians to be in range [-pi, pi].
 * @param angle in radians to be normalized.
 * @return normalized angle in range [-pi, pi].
 */
double Tools::normalizeAngle(double angle) {
    if (angle > M_PI || angle < -M_PI) {
        angle -= int(angle / M_PI) * M_PI;
    }
    return angle;
}
