#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

namespace Tools {

    /**
    * A helper method to calculate RMSE.
    */
    Eigen::VectorXd
    calculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

    /**
     * Normalizes angle in radians to be in range [-pi, pi].
     * @param angle in radians to be normalized.
     * @return normalized angle in range [-pi, pi].
     */
    double normalizeAngle(double angle);

    /**
     * Calculates NIS
     * @param S measurement covariance matrix
     * @param z_diff difference between predicted measurement and measurement vector
     * @return NIS value
     */
    double calculateNIS(const Eigen::MatrixXd& S, const Eigen::VectorXd& z_diff);
};

#endif /* TOOLS_H_ */
