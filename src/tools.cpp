#include <iostream>
#include <math.h>
#include "tools.h"


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

////////////////////////////////////////////////////////
// @function: calculate RMSE
// @parameters: vector of estimations, vector of ground_truths
// @return: rmse

// Evaluate the performance of a KF/EKF using rmse
// rmse = sqrt(sum((estimations - ground_truths)^2)/n)
// where n = size of estimations = size of ground_truth vectors
/////////////////////////////////////////////////////////
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// Ensure that the estimation vector size is non-zero
	if (estimations.size() == 0){
		cout << "Error! Size of estimations is zero!" << endl;
		return rmse;
	}

	// Ensure that the size of the estimations vector is the same as the size of the ground-truth vector
	if (estimations.size() != ground_truth.size()){
		cout << "Error! Size of estimations does not match size of ground truth!" << endl;
		return rmse;
	}

	// Loop through the estimations and accumulate residuals
	for(int i=0;i<estimations.size();++i){
		VectorXd residual = ground_truth[i] - estimations[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse/estimations.size();

	// Calculate the square root
	rmse = rmse.array().sqrt();
	cout << "RMSE = " << rmse << endl;
	return rmse;
}

////////////////////////////////////////////////////////
// @function: calculate h(x): the mapping from cartesian to polar coordinates
// @parameters: state vector
// @return: h(x) matrix
/////////////////////////////////////////////////////////
VectorXd Tools::CalculateHofX(const VectorXd& x_state){

	VectorXd H(3);

	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float rho = sqrt(px*px+py*py);
	float phi = atan2(py,px);
	float rho_dot;

	if(fabs(rho) < 0.0001)
		rho_dot = 0.0;
	else
		rho_dot = (px*vx+py*vy)/rho;

	H << rho, phi, rho_dot;

	return H;

}

////////////////////////////////////////////////////////
// @function: calculate Jacobian Hj for a nonlinear  h(x)
// @parameters: state vector
// @return: Jacobian matrix

/////////////////////////////////////////////////////////
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);

	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//Pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
