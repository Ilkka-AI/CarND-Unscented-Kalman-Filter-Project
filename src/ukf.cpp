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
  use_laser_ = true;//true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ =  true;

  ///* State dimension
  n_x_=5;

  ///* Augmented state dimension
  n_aug_=7;
 
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_,n_x_);
  
  is_initialized_= false;
  // Process noise standard deviation longitudinal acceleration in m/s^2

  std_a_ = 0.4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  
  std_yawdd_ = 0.2;

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
   
  time_us_=0;
  
  EPS= 0.002;

  // Counters for laser and radar measurements for collecting NIS values
  k_r=0;
  k_l=0;
   
  ///* Sigma point spreading parameter
  lambda_= 3 - n_aug_;
  

  Xsig_pred_= MatrixXd(n_x_, 2 * n_aug_ + 1);

  R_radar_ = MatrixXd(3, 3); 
  R_lidar_ = MatrixXd(2, 2);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;
 
  R_lidar_ << std_laspx_*std_laspx_,0,
0,std_laspy_*std_laspy_;

 
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  if(!is_initialized_){
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
     
    VectorXd weights = VectorXd(2*n_aug_+1);
    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_+lambda_);
      weights(i) = weight;
  }
   weights_=weights;
   
   x_ << 0,0,0	,0,0;
   P_<< 1,0,0,0,0,
       0,1,0,0,0,
       0,0,1,0,0,
       0,0,0,1,0,
       0,0,0,0,1; 
   if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
		x_(0) = EPS;
		x_(1) = EPS;
 }  
}
  
 else if(meas_package.sensor_type_ == MeasurementPackage::RADAR & use_radar_==true){
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    std::cout << "dt= " << std::endl << dt << std::endl;
    Prediction(dt);	
    UpdateRadar(meas_package);
    NIS_values_radar.resize(k_r+1);
    NIS_values_radar(k_r)=NIS;
    k_r++;
    //cout<<"NIS values radar" << NIS_values_radar << endl;
}
else if(meas_package.sensor_type_ == MeasurementPackage::LASER & use_laser_==true){
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    std::cout << "dt= " << std::endl << dt << std::endl;
    Prediction(dt);	
    UpdateLidar(meas_package);
    NIS_values_lidar.resize(k_l+1);
    NIS_values_lidar(k_l)=NIS;
    k_l++;

    //cout<<"NIS values lidar" << NIS_values_lidar << endl;
}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;  

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
 
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
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
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
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
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  x_=x;
  P_=P;
  //print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P << std::endl;
  Xsig_pred_=Xsig_pred;

}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
   
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
   
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  =  Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);   //phi                    
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
}

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
   S = S + R_radar_;

  
  //print result
  // std::cout << "Weights: " << std::endl << weights_ << std::endl;
  //std::cout << "Zsig: " << std::endl << Zsig << std::endl;
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  VectorXd z(3); 
  z << meas_package.raw_measurements_;

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  //std::cout << "z measurement: " << std::endl << z << std::endl;
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "z_diff: " << std::endl << z_diff << std::endl;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  NIS=z_diff.transpose()*Sinv*z_diff;


  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
 
 int n_z = 2;
   
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  =  Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = p_x;          
    Zsig(1,i) = p_y;   

}

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix 2x2
 
  S = S + R_lidar_;

  
  //print result
  //std::cout << "Weights: " << std::endl << weights_ << std::endl;
  //std::cout << "Zsig: " << std::endl << Zsig << std::endl;
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  VectorXd z(2); 

  z << meas_package.raw_measurements_;
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  //std::cout << "z measurement: " << std::endl << z << std::endl;
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "z_diff: " << std::endl << z_diff << std::endl;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  NIS=z_diff.transpose()*Sinv*z_diff;
  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}





