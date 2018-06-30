#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  //process noise covariance matrix.
  ekf_.Q_ = MatrixXd(4,4);

  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0,0,0,0;

  //state covariance Matrix . Q is 4x4 so is P since P'=F*P*F_t+Q
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  ekf_.F_ = MatrixXd(4,4);
  // dt is set to 1 for the time being.
  ekf_.F_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

 //Given .See prediction section of ProcessMeasurement();
  noise_ax = 9;
  noise_ay = 9;
  //since we are passing by address,the changes in Q_ should bre reflected in the values of ekf.Q_ ?
  // ekf_.Init(x_, P_, F_, H_, R_, Q_);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  float px=0,py=0,vx=0,vy=0;
  float dt=0.0;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      px = rho * cos(phi);
      py = rho * sin(phi);
      // Although radar gives velocity data in the form of the range rate ρ˙​,
      // a radar measurement does not contain enough information to determine the state variable velocities vx and vy
      vx = 0;
      vy = 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
      vx = 0;//no data ,so setting it to zero
      vy = 0;//no data ,so setting it to zero.

    }
    previous_timestamp_= measurement_pack.timestamp_;

    ekf_.x_<< px,py,vx,vy;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

   dt =(measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;
   ekf_.F_(0,2) = dt;
   ekf_.F_(1,3) = dt;
   ekf_.Q_<< pow(dt,4)*noise_ax/4,0,pow(dt,3)*noise_ax/2,0,
	         0,pow(dt,4)*noise_ay/4,0,pow(dt,3)*noise_ay/2,
	         pow(dt,3)*noise_ax/2,0,pow(dt,2)*noise_ax,0,
	         0,pow(dt,3)*noise_ay/2,0,pow(dt,2)*noise_ay;
  if (dt >0.001){
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
 if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
   ekf_.R_ = R_radar_;
   ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
   ekf_.UpdateEKF(measurement_pack.raw_measurements_);
 } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
   ekf_.R_ = R_laser_;
   ekf_.H_ = H_laser_;
   ekf_.Update(measurement_pack.raw_measurements_);
 }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
