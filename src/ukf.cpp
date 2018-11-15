#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
    std_a_ = 30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

    n_x_ = static_cast<int>(x_.size());

    n_aug_ = n_x_ + 2;

    lambda_ = 3 - n_x_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }
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

    //
    // 1. predict sigma points
    //

    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(n_x_) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // create augmented covariance matrix
    MatrixXd Q(2, 2);
    Q << std_a_ * std_a_, 0,
         0, std_yawdd_ * std_yawdd_;

    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug.bottomRightCorner(2, 2) = Q;

    // calculate the square root of the covariance matrix
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    double mag = sqrt(lambda_ + n_x_);
    Xsig_aug.col(0) = x_aug;
    for(int i = 0; i < n_aug_; ++i) {
        Xsig_aug.col(i + 1) = x_aug + mag * L.col(i);
        Xsig_aug.col(i + n_aug_ + 1) = x_aug - mag * L.col(i);
    }

    // create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {

        VectorXd col = Xsig_aug.col(i);

        double v = col(2);
        double yaw = col(3);
        double yaw_rate = col(4);
        double noise_a = col(5);
        double noise_yaw = col(6);

        VectorXd x = col.head(5);
        VectorXd m1 = VectorXd(5);
        if (fabs(yaw_rate) > 0.001) {
            double t1 = v / yaw_rate;
            double t2 = yaw + yaw_rate * delta_t;
            m1(0) = t1 * (sin(t2) - sin(yaw));
            m1(1) = t1 * (cos(yaw) - cos(t2));
        } else {
            m1(0) = v * cos(yaw) * delta_t;
            m1(1) = v * sin(yaw) * delta_t;
        }
        m1(2) = 0;
        m1(3) = yaw_rate * delta_t;
        m1(4) = 0;

        VectorXd m2 = VectorXd(5);
        double t3 = 0.5 * (delta_t*delta_t);
        m2 << t3 * cos(yaw) * noise_a,
                t3 * sin(yaw) * noise_a,
                delta_t * noise_a,
                t3 * noise_yaw,
                delta_t * noise_yaw;

        Xsig_pred.col(i) = x + m1 + m2;
    }

    // predict mean and covariance

    VectorXd weights = VectorXd(2 * n_aug_ + 1);

    // set weights
    weights(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
        weights(i) = 1.0 / (2 * (lambda_ + n_aug_));
    }

    //
    // 2. predict state mean
    //
    x_.fill(0.0);
    for(int i = 0; i < n_x_; ++i) {
        x_ += weights(i) * Xsig_pred.col(i);
    }

    //
    // 3. predict state covariance matrix
    //
    P_.fill(0.0);
    for(int i = 0; i < 2 * n_aug_ + 1; ++i) {
        MatrixXd col = Xsig_pred.col(i) - x_;

        // normalizing the angle
        while (col(3)> M_PI) col(3)-=2.*M_PI;
        while (col(3)<-M_PI) col(3)+=2.*M_PI;

        MatrixXd c2 = col.transpose();
        P_ += weights(i) * col * c2;
    }
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
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
}
