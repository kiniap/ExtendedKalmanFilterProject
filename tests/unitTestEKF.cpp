/*
 * rmseTest.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: kiniap
 */

#define BOOST_TEST_MODULE unitTestsEKF
#include <boost/test/unit_test.hpp>
#include "Eigen/Dense"
#include "../src/tools.h"
#include <fstream>
#include <iostream>
#include "../src/measurement_package.h"
#include "tracking.h"
#include "../src/FusionEKF.h"

using namespace std;

const double tolerance = 0.00001;

/*
  Test the computation of the RMSE
 */
BOOST_AUTO_TEST_CASE(RMSE_TEST)
{
	vector<VectorXd> estimations;
	vector<VectorXd> ground_truth;

	//the input list of estimations
	VectorXd e(4);
	e << 1, 1, 0.2, 0.1;
	estimations.push_back(e);
	e << 2, 2, 0.3, 0.2;
	estimations.push_back(e);
	e << 3, 3, 0.4, 0.3;
	estimations.push_back(e);

	//the corresponding list of ground truth values
	VectorXd g(4);
	g << 1.1, 1.1, 0.3, 0.2;
	ground_truth.push_back(g);
	g << 2.1, 2.1, 0.4, 0.3;
	ground_truth.push_back(g);
	g << 3.1, 3.1, 0.5, 0.4;
	ground_truth.push_back(g);

	Tools ekf_tools;
	//call the CalculateRMSE
	VectorXd rmse = ekf_tools.CalculateRMSE(estimations, ground_truth);
	VectorXd rmseResult(4);
	rmseResult << 0.1,0.1,0.1,0.1;
	bool r = (rmse-rmseResult).norm() < tolerance;
	BOOST_CHECK_EQUAL(true, r);
}


/*
 Test the Computation of the Jacobian Matrix
*/
BOOST_AUTO_TEST_CASE(JACOBIAN_TEST)
{
	//predicted state  example
	//px = 1, py = 2, vx = 0.2, vy = 0.4
	VectorXd x_predicted(4);
	x_predicted << 1, 2, 0.2, 0.4;

	// Call the calcuateJacobian functions
	Tools ekf_tools;
	MatrixXd Hj = ekf_tools.CalculateJacobian(x_predicted);
	MatrixXd HjResult(3,4);
	HjResult << 0.447214, 0.894427, 0, 0,
			    -0.4, 0.2, 0, 0,
				0, 0, 0.447214, 0.894427;

	bool r = (Hj-HjResult).norm() < tolerance;
	BOOST_CHECK_EQUAL(true, r);
}


BOOST_AUTO_TEST_CASE(KALMAN_FILTER)
{
	vector<MeasurementPackage> measurement_pack_list;

	// hardcoded input file with laser and radar measurements
	string in_file_name_ = "obj_pose-laser-radar-synthetic-input.txt";
	ifstream in_file(in_file_name_.c_str(),std::ifstream::in);

	if (!in_file.is_open()) {
		cout << "Cannot open input file: " << in_file_name_ << endl;
	}

	string line;
	// set i to get only first 3 measurments
	int i = 0;
	while(getline(in_file, line) && (i<=3)){

		MeasurementPackage meas_package;

		istringstream iss(line);
		string sensor_type;
		iss >> sensor_type;	//reads first element from the current line
		int64_t timestamp;
		if(sensor_type.compare("L") == 0){	//laser measurement
			//read measurements
			meas_package.sensor_type_ = MeasurementPackage::LASER;
			meas_package.raw_measurements_ = VectorXd(2);
			float x;
			float y;
			iss >> x;
			iss >> y;
			meas_package.raw_measurements_ << x,y;
			iss >> timestamp;
			meas_package.timestamp_ = timestamp;
			measurement_pack_list.push_back(meas_package);

		}else if(sensor_type.compare("R") == 0){
			//Skip Radar measurements
			continue;
		}
		i++;

	}

	//Create a Tracking instance
	Tracking tracking;
	VectorXd x;
	//call the ProcessingMeasurement() function for each measurement
	size_t N = measurement_pack_list.size();
	for (size_t k = 0; k < N; ++k) {	//start filtering from the second frame (the speed is unknown in the first frame)
		x = tracking.ProcessMeasurement(measurement_pack_list[k]);

	}

	if(in_file.is_open()){
		in_file.close();
	}

	VectorXd x_final(4);
	x_final << 1.34291, 0.364408,  2.32002, -0.722813;


	bool r = (x - x_final).norm() < tolerance;
	BOOST_CHECK_EQUAL(true, r);

}

BOOST_AUTO_TEST_CASE(EKF)
{
	//vector<MeasurementPackage> measurement_pack_list;
	Tools tools;
	vector<VectorXd> estimations;
	vector<VectorXd> ground_truth;

	VectorXd RMSE;

	//Create a FusionEKF Insstance
	FusionEKF fusionEKF;

	// hardcoded input file with laser and radar measurements
	string in_file_name_ = "obj_pose-laser-radar-synthetic-input.txt";
	ifstream in_file(in_file_name_.c_str(),std::ifstream::in);

	if (!in_file.is_open()) {
		cout << "Cannot open input file: " << in_file_name_ << endl;
	}

	string line;
	while(getline(in_file, line)){

		MeasurementPackage meas_package;

		istringstream iss(line);
		string sensor_type;
		iss >> sensor_type;	//reads first element from the current line
		int64_t timestamp;
		if(sensor_type.compare("L") == 0){	//laser measurement
			//continue; // uncomment this line to ignore laser measurements
			meas_package.sensor_type_ = MeasurementPackage::LASER;
			meas_package.raw_measurements_ = VectorXd(2);
			float x;
			float y;
			iss >> x;
			iss >> y;
			meas_package.raw_measurements_ << x,y;
			iss >> timestamp;
			meas_package.timestamp_ = timestamp;
			//measurement_pack_list.push_back(meas_package);

		}else if (sensor_type.compare("R") == 0) {
			//continue; // uncomment this line to ignore radar measurements
  	  		meas_package.sensor_type_ = MeasurementPackage::RADAR;
      		meas_package.raw_measurements_ = VectorXd(3);
      		float ro;
  	  		float theta;
  	  		float ro_dot;
      		iss >> ro;
      		iss >> theta;
      		iss >> ro_dot;
      		meas_package.raw_measurements_ << ro,theta, ro_dot;
      		iss >> timestamp;
      		meas_package.timestamp_ = timestamp;
			//measurement_pack_list.push_back(meas_package);
      }

      float x_gt;
  	  float y_gt;
  	  float vx_gt;
  	  float vy_gt;
  	  iss >> x_gt;
  	  iss >> y_gt;
  	  iss >> vx_gt;
  	  iss >> vy_gt;
  	  VectorXd gt_values(4);
  	  gt_values(0) = x_gt;
  	  gt_values(1) = y_gt;
  	  gt_values(2) = vx_gt;
  	  gt_values(3) = vy_gt;
  	  ground_truth.push_back(gt_values);

      //Call ProcessMeasurment(meas_package) for Kalman filter
	  fusionEKF.ProcessMeasurement(meas_package);

	  //Push the current estimated x,y positon from the Kalman filter's state vector

	  VectorXd estimate(4);

	  double p_x = fusionEKF.ekf_.x_(0);
	  double p_y = fusionEKF.ekf_.x_(1);
	  double v1  = fusionEKF.ekf_.x_(2);
	  double v2 = fusionEKF.ekf_.x_(3);

	  estimate(0) = p_x;
	  estimate(1) = p_y;
	  estimate(2) = v1;
	  estimate(3) = v2;

	  estimations.push_back(estimate);

	  RMSE = tools.CalculateRMSE(estimations, ground_truth);
	}


	if(in_file.is_open()){
		in_file.close();
	}

	//px, py, vx, vy output coordinates must have an RMSE <= [.11, .11, 0.52, 0.52]
    //when using the file: "obj_pose-laser-radar-synthetic-input.txt which is the same data file the simulator uses for Dataset 1"
	bool r = (RMSE(0) < 0.11) && (RMSE(1) < 0.11) && (RMSE(2) < 0.52) && (RMSE(3) < 0.52);
	BOOST_CHECK_EQUAL(true, r);

}

