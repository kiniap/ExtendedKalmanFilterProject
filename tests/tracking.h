
#ifndef TRACKING_H_
#define TRACKING_H_

#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include "../src/measurement_package.h"
#include "../src/kalman_filter.h"

class Tracking {
public:
	Tracking();
	virtual ~Tracking();
	Eigen::VectorXd ProcessMeasurement(const MeasurementPackage &measurement_pack);
	KalmanFilter kf_;

private:
	bool is_initialized_;
	int64_t previous_timestamp_;

	//acceleration noise components
	float noise_ax;
	float noise_ay;

};

#endif /* TRACKING_H_ */
