/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (is_initialized){
		return;
	}

	num_particles=10;

	// For readability
	double std_x=std[0];
	double std_y=std[1];
	double std_theta=std[2];

	// Create normal distribution for x,y,theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	default_random_engine gen;
	//Initliaze particles
	for (int i = 0; i < num_particles; ++i) {
		particles.push_back(Particle());
		particles[i].id=i;
		particles[i].x=dist_x(gen);
		particles[i].y=dist_y(gen);
		particles[i].theta=dist_theta(gen);
		particles[i].weight=1;
		//create a weight for each particle initialized
		weights.push_back(1.0);
	}
	 is_initialized = true;
	 cout<<"initialized";
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// For readability
	double std_x=std_pos[0];
	double std_y=std_pos[1];
	double std_theta=std_pos[2];

	// Create normal distribution for x,y,theta with 0 mean.
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);
	default_random_engine gen;
	for ( int i=0;i<num_particles;i++){
		if (fabs(yaw_rate)!=0){
			particles[i].x+=(velocity/yaw_rate)*(sin(particles[i].theta+ (yaw_rate*delta_t))-sin(particles[i].theta));
			particles[i].y+=(velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta+(yaw_rate*delta_t)));
			particles[i].theta+=yaw_rate*delta_t;
		}else{
			// very little or no change in yaw(theta)
			particles[i].x+=velocity*cos(particles[i].theta)*delta_t;
			particles[i].y+=velocity*sin(particles[i].theta)*delta_t;
			// yaw(theta) remains the same
		}


		//Add random noise
		particles[i].x+=dist_x(gen);
		particles[i].y+=dist_y(gen);
		particles[i].theta+=dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	// Iterate over observations.These might be less than equal to the number of predicted measurements
	for (int i=0;i<observations.size();i++){
		double min_dist_sq=numeric_limits<double>::max();
		int min_id=-1;
		for (int j=0;j<predicted.size();j++){
			double x_diff=predicted[j].x-observations[i].x;
			double y_diff=predicted[j].y-observations[i].y;
			double dist_sq=(x_diff*x_diff) +(y_diff*y_diff);
			if (dist_sq < min_dist_sq){
				min_dist_sq=dist_sq;
				min_id=predicted[j].id;
			}
		}
		observations[i].id=min_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	int num_obs=observations.size();
	int num_landmarks=map_landmarks.landmark_list.size();
	std::vector<LandmarkObs> observations_map_coords(num_obs);
	std::vector<LandmarkObs> landmarks_inrange;



	//Iterate over all particles
	for (int i = 0;i < num_particles;i++){
		particles[i].weight = 1.0;


		//Iterate over all landmarks to see whats in range of this particle
		for (int k = 0;k < num_landmarks;k++){
			double x_diff = map_landmarks.landmark_list[k].x_f - particles[i].x;
			double y_diff = map_landmarks.landmark_list[k].y_f - particles[i].y;

			if (((x_diff*x_diff)+(y_diff*y_diff)) <= sensor_range*sensor_range){
				LandmarkObs landmark;
				landmark.id=map_landmarks.landmark_list[k].id_i;
				landmark.x=map_landmarks.landmark_list[k].x_f;
				landmark.y=map_landmarks.landmark_list[k].y_f;
				landmarks_inrange.push_back(landmark);
			}
		}
		double x_p = particles[i].x;
		double y_p = particles[i].y;
		double theta = particles[i].theta;

		// Convert all observations to MAP coordinates for this particle.
		for (int j=0;j<num_obs;j++){
			//car observation coordinates x_c ​ and y_c ,map particle coordinates x_p ​ and y_p, and
			// rotation angle (theta)

			double x_c = observations[j].x;
			double y_c = observations[j].y;
			observations_map_coords[j].id = observations[j].id;
			observations_map_coords[j].x = x_p + cos(theta)*x_c - sin(theta)*y_c;
			observations_map_coords[j].y = y_p + sin(theta)*x_c + cos(theta)*y_c;
		}

		// Calculate predicted distances to in-range landmarks.
		int num_landmarks_inrange=landmarks_inrange.size();

		// Populate the nearest landmark ID in observations_map_coords
		dataAssociation(landmarks_inrange,observations_map_coords);

		//Calculate weights for the particle.

		//for each observation of that particle ,see which landmark it is mapped to and update weight
		for (int j=0;j<num_obs;j++){
			double x_diff,y_diff;
			// Iterate through all in range landmarks to see which one this observation is mapped to.
			for (int k=0;k<num_landmarks_inrange;k++){
				if (observations_map_coords[j].id==landmarks_inrange[k].id)
				{
					x_diff=observations_map_coords[j].x-landmarks_inrange[k].x;
					y_diff=observations_map_coords[j].y-landmarks_inrange[k].y;
					break;// once you find a match break
				}
			}
			//Absolute of the covariance matrix = (variance_xx*variance_yy - 0);
			//variance_xx =std_x*std_x
			double std_x=std_landmark[0];
			double std_y=std_landmark[1];
			double den=sqrt(2*M_PI)*std_x*std_y;
			double num=exp(-0.5*((x_diff*x_diff/(std_x*std_x))+((y_diff*y_diff)/(std_y*std_y))));
			double weight = num/den;
			//Don't want to make particle probability 0.
			if (weight==0){
				weight=0.00001;
			}
			particles[i].weight*=weight;
		}
		// Update vector that keeps tracks of weights of all particles.It's used in resampling .
		weights[i]=particles[i].weight;

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;

	// produces random integers on the interval [0, n), where the probability of each individual integer i is defined as
	//wi/S, that is the weight of the ith integer divided by the sum of all n weights.
	//http://www.cplusplus.com/reference/random/discrete_distribution/discrete_distribution/
  std::discrete_distribution<> d (weights.begin(),weights.end());
  // initialise new particle array
  vector<Particle> resampledParticles;
  // resample particles
  for (int i = 0; i < num_particles; i++){
		resampledParticles.push_back(particles[d(gen)]);
	}
	particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
