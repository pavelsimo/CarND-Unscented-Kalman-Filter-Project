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
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    // if this is false, lidar measurements will be ignored (except during init)
    use_lidar_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.231;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.20;

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

    // set weights
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
        double weight = 0.5 / (lambda_ + n_aug_);
        weights_(i) = weight;
    }

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    // custom variables

    use_nis_ = true;
    nis_ = 0;
    if (use_nis_) {
        cerr << "sensor_type" << "," << "timestamp" << "," << "nis" << '\n';
    }
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
    cout << "Start - ProcessMeasurement" << endl;
    if (!is_initialized_) {

        VectorXd z = meas_package.raw_measurements_;

        P_.fill(0.0);
        P_(0, 0) = 0.30; // px
        P_(1, 1) = 0.30; // py
        P_(2, 2) = 1; // v (magnitude)
        P_(3, 3) = 1; // yaw
        P_(4, 4) = 1; // yaw rate

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            // convert radar from polar to cartesian coordinates and initialize state.
            // radial distance from the origin
            double ro = z[0];

            // angle between rho and x
            double theta = z[1];

            x_ << ro * cos(theta), ro * sin(theta), 0, 0, 0;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            double px = z[0];
            double py = z[1];
            x_ << px, py, 0, 0, 0;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;
        cout << "Done - ProcessMeasurement" << endl;
        return;
    }

    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    Prediction(delta_t);

    char sensor_type = 'X';
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        if (use_lidar_) {
            UpdateLidar(meas_package);
            sensor_type = 'L';
        }
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        if (use_radar_) {
            UpdateRadar(meas_package);
            sensor_type = 'R';
        }
    }
    time_us_ = meas_package.timestamp_;
    if (use_nis_) {
        cerr << sensor_type << "," << time_us_ << "," << nis_ << '\n';
    }
    cout << "x = " << x_ << endl;
    cout << "P = " << P_ << endl;
    cout << "Done - ProcessMeasurement" << endl;
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
    cout << "Start - Prediction" << endl;
    cout << "delta_t = " << delta_t << endl;

    // ===============================================
    // 1. create augmented sigma points
    // ===============================================

    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    // create augmented mean state
    x_aug.head(n_x_) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    // create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    // ===============================================
    // 2. predict sigma points
    // ===============================================

    // predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        double px = Xsig_aug(0, i);
        double py = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawrate = Xsig_aug(4, i);
        double noise_a = Xsig_aug(5, i);
        double noise_yawrate = Xsig_aug(6, i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawrate) > 0.001) {
            px_p = px + v / yawrate * (sin(yaw + yawrate * delta_t) - sin(yaw));
            py_p = py + v / yawrate * (cos(yaw) - cos(yaw + yawrate * delta_t));
        } else {
            px_p = px + v * delta_t * cos(yaw);
            py_p = py + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawrate * delta_t;
        double yawrate_p = yawrate;

        // add noise
        px_p = px_p + 0.5 * noise_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * noise_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + noise_a * delta_t;

        yaw_p = yaw_p + 0.5 * noise_yawrate * delta_t * delta_t;
        yawrate_p = yawrate_p + noise_yawrate * delta_t;

        // write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawrate_p;
    }

    // ===============================================
    // 3. predict state mean and covariance
    // ===============================================

    // predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    // predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        Tools::NormalizeAngle(x_diff(3));

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }

    cout << "Done - Prediction" << endl;
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
    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_pred = x_.head(2);
    MatrixXd H(2, n_x_);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;
    MatrixXd R(2, 2);
    R << std_laspx_ * std_laspx_, 0,
         0, std_laspy_ * std_laspy_;

    VectorXd y = z - H * x_;
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
    P_ = (I - K * H) * P_;

    nis_ = Tools::NIS(S, z_pred, z);
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

    cout << "Start - UpdateRadar" << endl;

    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        // r
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);
        // phi
        Zsig(1, i) = atan2(p_y, p_x);
        // r_dot
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);
    }

    // calculate mean predicted measurement
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        z_pred += weights_(i) * Zsig.col(i);
        Tools::NormalizeAngle(z_pred(1));
    }

    // calculate innovation covariance matrix S

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        MatrixXd z_diff = Zsig.col(i) - z_pred;
        Tools::NormalizeAngle(z_diff(1));

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
    S = S + R;

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        Tools::NormalizeAngle(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        Tools::NormalizeAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K
    MatrixXd K = Tc * S.inverse();

    // update state mean and covariance matrix
    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_diff = z - z_pred;
    Tools::NormalizeAngle(z_diff(1));

    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
    nis_ = Tools::NIS(S, z_pred, z);

    cout << "Done - UpdateRadar" << endl;
}
