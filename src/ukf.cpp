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
    
  // Print NIS
    
    print_nis_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2 0,2
    std_a_ = 0.6;

  // Process noise standard deviation yaw acceleration in rad/s^2 0,19
    std_yawdd_ = 0.6; //1.2;  // 0.5

  // Laser measurement noise standard deviation position1 in m  0.15
    std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m 0.3
    std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad 0.03
    std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    is_initialized_ = false;
    n_x_ = 5;
    n_aug_ = 7;
    n_z_ = 3;
    n_l_ = 2;
    lambda_ = 3 - n_aug_;
    weights_ = Eigen::VectorXd(2*n_aug_+1);
    
    /* Initializa weights */
    weights_.fill(1.0 / 2.0 / (lambda_ + n_aug_));
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    
    R_radar_ = MatrixXd(n_z_, n_z_);
    R_radar_.fill(0);
    R_radar_(0,0) = std_radr_*std_radr_;
    R_radar_(1,1) = std_radphi_*std_radphi_;
    R_radar_(2,2) = std_radrd_*std_radrd_;

    R_laser_ = MatrixXd(2, n_l_);
    R_laser_.fill(0);
    R_laser_(0,0) = std_laspx_*std_laspx_;
    R_laser_(1,1) = std_laspy_*std_laspy_;
 
    
    
    
    H_ = MatrixXd(2, 5);
    H_.fill(0);
    H_(0, 0) = 1;
    H_(1,1) = 1;

}

UKF::~UKF() {}

double UKF::Normalize(double value){
    return atan2(sin(value), cos(value));

}

void UKF::Initialize(MeasurementPackage measurement_pack){
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement

        x_ = VectorXd(5);
        x_ << 1, 1, 0, 0, 0;
    
        time_us_ = measurement_pack.timestamp_;
        counter_ = 0;
        
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            
            double r = measurement_pack.raw_measurements_[0];
            double theta = measurement_pack.raw_measurements_[1];
            
            x_(0)  = r * cos(theta);
            x_(1)  = r * sin(theta);
            
            P_ = MatrixXd(5,5);
            
            P_ << pow(R_radar_(0,0)*cos(theta)+r*sin(theta)*R_radar_(1,1),2), 0, 0, 0, 0,
            0, pow(R_radar_(0,0)*sin(theta)+r*cos(theta)*R_radar_(1,1),2), 0, 0,0,
            //P_ << 1, 0, 0, 0, 0,
            //0, 1, 0, 0,0,
            0, 0, 1, 0,0,
            0, 0, 0, 1,0,
            0, 0, 0, 0, 1;
            
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            
            x_[0] = measurement_pack.raw_measurements_[0];
            x_[1] = measurement_pack.raw_measurements_[1];
            x_[2] = 0;
            x_[3] = 0; 
            x_[4] = 0;
            
            P_ = MatrixXd(5,5);
            
            P_ << std_laspx_*std_laspx_, 0, 0, 0,0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            //P_ << 1, 0, 0, 0,0,
            //0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
            
        }
        
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    
    if (!is_initialized_){
        Initialize(meas_package);
        return;
    }

    double delta_t = ( meas_package.timestamp_ - time_us_) / 1000000.0;
    double nis;
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
        counter_ += 1;
        std::tie(x_, P_) = Prediction(delta_t);
        std::tie(x_, P_, nis) = UpdateLidar(meas_package, x_, P_, R_laser_, H_, n_x_);
        if (print_nis_){
            cout << "Laser_NIS \t" << nis << "\t" << x_(3) << endl;
        }

        //cout << "Updated_Lidar_x_" << x_.transpose() << endl;
        //cout << "Updated Liudar P_" << P_ << endl;

        time_us_ = meas_package.timestamp_;
        
    }else if (meas_package.sensor_type_== MeasurementPackage::RADAR && use_radar_){
        counter_ += 1;
        std::tie(x_, P_) = Prediction(delta_t);
        //UpdateRadar(meas_package);
        std::tie(x_, P_, nis) = UpdateRadar(meas_package, x_, P_, Xsig_pred_, weights_, R_radar_, n_x_, n_z_, n_aug_);
        if (print_nis_){
            cout << "Radar_NIS \t" << nis << "\t" << x_(3)<< endl;
        }
        //cout << "Updated_Radar_x_" << x_.transpose() << endl;
        //cout << "Updated Radar P_" << P_ << endl;
        
        time_us_ = meas_package.timestamp_;
    }
 }
/**
 * ComputeSigmaPoints takes actual state and covariance matrix
 * and computes sigma points
 */
MatrixXd UKF::ComputeSigmaPoints(VectorXd x, MatrixXd P, int n_x, int n_aug, double lambda, double std_a, double std_yawdd){
    // Define augmented variables
    
    VectorXd x_aug = VectorXd(n_aug);
    MatrixXd P_aug = MatrixXd(n_aug, n_aug);
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2*n_aug+1);
    
    
    // Compute augmented values
    x_aug.fill(0);
    x_aug.head(n_x) = x;
    
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x, n_x) = P;
    P_aug(n_x,n_x) = std_a*std_a;
    P_aug(n_x+1,n_x+1) = std_yawdd*std_yawdd;
    
    // Compute sigma points
    MatrixXd L = P_aug.llt().matrixL();
    
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug; i++)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
        //Xsig_aug(3,i+1) = Normalize(Xsig_aug(3,i+1));
        //Xsig_aug(3,i+1+n_aug_) = Normalize(Xsig_aug(3,i+1+n_aug_));
    }
    return Xsig_aug;
}

/**
 * PredictSigmaPoints computes new predicted Sigma Points
 * and returns them
 */

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t, int n_x, int n_aug){
    
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    
    for (int i = 0; i< 2*n_aug+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        //yaw_p = Normalize(yaw_p);
        yawd_p = yawd_p + nu_yawdd*delta_t;
 
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
        
    }
    
    return Xsig_pred;
    
}

/**
 * ComputePredictedState computes new predicted state and covariance matrix
 * from  predicted sigma points
 */

std::tuple<VectorXd, MatrixXd> UKF::ComputePredictedState(MatrixXd Xsig_pred, VectorXd weights, int n_x, int n_aug){
    
    VectorXd x = VectorXd(n_x);
    x.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        x = x + weights(i) * Xsig_pred.col(i);
    }
    x(3) = Normalize(x(3)); // Just in case

    //predicted state covariance matrix
    MatrixXd P = MatrixXd(n_x, n_x);
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        VectorXd x_diff = Xsig_pred.col(i) - x;
        x_diff(3) = Normalize(x_diff(3));
        P = P + weights(i) * x_diff * x_diff.transpose() ;
    }
    

    return std::make_tuple(x, P);
}


/**
 * ComputePredictedMeasurementState computes new predicted state and covariance matrix
 * from  predicted sigma points for Radar measure
 * vectors have r, phi, rdot
 */

std::tuple<VectorXd, MatrixXd> UKF::ComputePredictedMeasurementdStateRadar(MatrixXd Zsig, VectorXd weights, MatrixXd R, int n_z, int n_aug){
    
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    
    // Compute mean
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    z_pred(1) = Normalize(z_pred(1));
    //predicted state covariance matrix
    
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        // state difference
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = Normalize(z_diff(1));
        S = S + weights(i) * z_diff * z_diff.transpose() ;
    }
    
    S = S + R;

    return std::make_tuple(z_pred, S);
    
}

/**
 * TransformSigmaPoints converts sigma points to radar measurement space Z
 * This is the Radar Measurement Model
 * Radar gives r, phi and rdot
 */

Eigen::MatrixXd  UKF::TransformSigmaPointsRadar(MatrixXd Xsig, int n_z, int n_aug){

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig(0,i);
        double p_y = Xsig(1,i);
        double v  = Xsig(2,i);
        double yaw = Xsig(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r

        //TODO: If p_x and p_y == 0 we need to get different computation

        Zsig(1,i) = atan2(p_y,p_x);
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    return Zsig;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
std::tuple<VectorXd, MatrixXd> UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    // Now we predict new state
    // First compute predicted sigma points
    //cout << "x_ before prediction" << endl << x_ << endl;
    //cout << "P_ before prediction" << endl << P_ << endl;
  
    MatrixXd Xsig_aug = ComputeSigmaPoints(x_, P_, n_x_, n_aug_, lambda_, std_a_, std_yawdd_);
    Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t, n_x_, n_aug_);
    return ComputePredictedState(Xsig_pred_, weights_, n_x_, n_aug_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
std::tuple<VectorXd, MatrixXd, double> UKF::UpdateLidar(MeasurementPackage meas_package, VectorXd x, MatrixXd P, MatrixXd R_laser, MatrixXd H, int n_x) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
   
   */
 
    VectorXd z = meas_package.raw_measurements_;
    VectorXd error = z - H * x;   // Compute error
 
    MatrixXd I = MatrixXd::Identity(n_x, n_x);
    
    MatrixXd Ht = H.transpose();
    MatrixXd S = (H * P * Ht) + R_laser;
    
    //cout << " P= " << P << endl;
    //cout << " S= " << S << endl;
    
    MatrixXd K = P * Ht * S.inverse();
    
    // Compute the nIS error
    
    double nis = error.transpose() * S.inverse() * error;
    //cout << " K= " << K << endl;
    
    VectorXd delta = K * error;
    x = x + delta;
    
    x(3) = Normalize(x(3));
    P = (I - K * H) * P;
    return std::make_tuple(x, P, nis);
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
std::tuple<VectorXd, MatrixXd, double> UKF::UpdateRadar(MeasurementPackage meas_package, VectorXd x,
                                                MatrixXd P, MatrixXd Xsig_pred, VectorXd weights,
                                                MatrixXd R_radar, int n_x, int n_z, int n_aug) {
    
    VectorXd z_pred;
    MatrixXd S;
    
    VectorXd z = meas_package.raw_measurements_;
    MatrixXd Zsig = TransformSigmaPointsRadar(Xsig_pred, n_z, n_aug);
    std::tie(z_pred, S) = ComputePredictedMeasurementdStateRadar(Zsig, weights, R_radar, n_z, n_aug);
    
    // Compute Tc Correlation Matrix
    
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = Normalize(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        x_diff(3) = Normalize(x_diff(3));
        
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }
    
    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;
    z_diff(1) = Normalize(z_diff(1));
    
    double nis = z_diff.transpose() * S.inverse() * z_diff;
    
    //update state mean and covariance matrix
    
    VectorXd delta = K * z_diff;
    x = x + delta;
    P = P - K*S*K.transpose();
    x(3) = Normalize(x(3));

    if (fabs(x(3)) > M_PI){
        cout << "Yaw out of margins when updating Radar" << x(3) << endl;
    }

    return std::make_tuple(x, P, nis);
}



