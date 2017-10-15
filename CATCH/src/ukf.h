#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

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
    
    bool print_nis_;
    
    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;
    
    ///* time when the state is true, in us
    long long time_us_;
    
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
    double std_radrd_ ;
    
    ///* Weights of sigma points
    VectorXd weights_;
    
    ///* State dimension
    int n_x_;
    
    ///* Radar measurement dimension
    int n_z_;
    
    ///* Laser measurement dimension
    
    int n_l_;
    
    ///* R is the Radar masurement covariances
    
    MatrixXd R_radar_;
    
    ///* R is the laser masurement covariances
    
    MatrixXd R_laser_;
    
    
    ///* H Laser
    
    MatrixXd H_;
    
    
    ///* Augmented state dimension
    int n_aug_;
    
    ///* Sigma point spreading parameter
    double lambda_;
    
    ///* Iteration counter
    
    int counter_;
    
    /**
     * Constructor
     */
    UKF();
    
    /**
     * Destructor
     */
    virtual ~UKF();
    
    /** Normalize
     *
     *  Normalizes angle variable
     */
    
    double Normalize(double value);
    /**
     * Initialize variables with first measurement
     */
    
    void Initialize(MeasurementPackage meas_package);
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
    std::tuple<VectorXd, MatrixXd> Prediction(double delta_t);
    
    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    std::tuple<VectorXd, MatrixXd, double> UpdateLidar(MeasurementPackage meas_package, VectorXd x, MatrixXd P, MatrixXd R_laser,  MatrixXd H,  int n_x );
    
    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    std::tuple<VectorXd, MatrixXd, double> UpdateRadar(MeasurementPackage meas_package, VectorXd x, MatrixXd P, MatrixXd Xsig_pred, VectorXd weights, MatrixXd R_radar, int n_x, int n_z, int n_aug);
    
    
    /**
     * ComputeSigmaPoints computes the Sigma Points matrix
     */
    
    MatrixXd ComputeSigmaPoints(VectorXd x, MatrixXd P, int n_x, int n_aug, double lambda, double std_a, double std_yawdd);
    
    /**
     * PredictSigmaPoints computes new predicted Sigma Points
     * and returns them
     */
    
    Eigen::MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t, int n_x, int n_aug);
    
    /**
     * ComputePredictedState computes new predicted state and covariance matrix
     * from  predicted sigma points
     */
    
    std::tuple<VectorXd, MatrixXd> ComputePredictedState(MatrixXd Xsig_pred, VectorXd weights, int n_x, int n_aug);
    
    /**
     * TransformSigmaPoints converts sigma points to radar measurement space Z
     */
    Eigen::MatrixXd  TransformSigmaPointsRadar(MatrixXd Xsig, int n_z, int n_aug);
    Eigen::MatrixXd  TransformSigmaPointsLidar(MatrixXd Xsig);
    /**
     * ComputePredictedMeasurementState computes new predicted state and covariance matrix
     * from  predicted sigma points
     * We could use the other function but angle normalization breaks it
     */
    
    std::tuple<VectorXd, MatrixXd> ComputePredictedMeasurementdStateRadar(MatrixXd Zsig, VectorXd weights, MatrixXd R,  int n_z, int n_aug);
    
};
#endif /* UKF_H */

