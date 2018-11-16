#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    /**
     * A helper method for angle normalization
     */
    static void NormalizeAngle(double &angle);

    /**
    * A helper method to calculate Normalized Innovation Squared (NIS)
    */
    static double NIS(const MatrixXd &S, const VectorXd &z_pred, const VectorXd &z);
};

#endif /* TOOLS_H_ */