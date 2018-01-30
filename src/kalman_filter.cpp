#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	x_ = F_*x_;
	P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	// Measurement update
	VectorXd y = z - H_*x_;

	UpdateCommon(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	Tools ekfTools;
	VectorXd h_of_x(3);
	h_of_x = ekfTools.CalculateHofX(x_);

	// Measurement update
	VectorXd y = z - h_of_x;
	NormalizeAngle(y(1));

	UpdateCommon(y);
}

void KalmanFilter::UpdateCommon(const VectorXd& y)
{
  const MatrixXd PHt = P_ * H_.transpose();
  const MatrixXd S = H_ * PHt + R_;
  const MatrixXd K = PHt * S.inverse();

  x_ += K * y;
  P_ -= K * H_ * P_;
}

void KalmanFilter::NormalizeAngle(double& phi)
{
	phi = atan2(sin(phi), cos(phi));
}
