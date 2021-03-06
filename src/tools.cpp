#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        return rmse;
    }

    for (int i = 0; i < estimations.size(); ++i) {
        VectorXd c1 = estimations[i] - ground_truth[i];
        c1 = c1.array() * c1.array();
        rmse += c1;
    }

    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

void Tools::NormalizeAngle(double &angle) {
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
}

double Tools::NIS(const MatrixXd &S, const VectorXd &z_pred, const VectorXd &z) {
    VectorXd z_diff = z - z_pred;
    double res = z_diff.transpose() * S.inverse() * z_diff;
    return res;
}