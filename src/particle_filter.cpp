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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	if (this->initialized() != true) //if didn't initialize
	{

		this->num_particles = 100;//how many particles do you want, I tried 50,100,200 and all of them work well.

		normal_distribution<double> dist_x(x, std[0]);
		normal_distribution<double> dist_y(y, std[1]);
		normal_distribution<double> dist_theta(theta, std[2]);
		default_random_engine gen; 

		Particle p;
		for (int i = 0; i < ParticleFilter::num_particles; i++)
		{
			p.id = i;
			p.x = dist_x(gen);
			p.y = dist_y(gen);
			p.theta = dist_theta(gen);
			//Note : don't use  p.x += dist_x(gen)
			p.weight = 1.0;
			this->particles.push_back(p);
			this->weights.push_back(p.weight);
		}
		this->is_initialized = true;
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	default_random_engine gen;
	for (int i = 0; i < this->num_particles; i++)
	{
		double x_old = this->particles[i].x;
		double y_old = this->particles[i].y;
		double yaw_old = this->particles[i].theta;
		double x_new;
		double y_new;
		double yaw_new;
		// check if yaw_rate == 0 or not, apply two different situations
		if (yaw_rate != 0)
		{
			x_new = x_old + velocity / (yaw_rate)*(sin(yaw_old + yaw_rate*delta_t) - sin(yaw_old));
			y_new = y_old + velocity / (yaw_rate)*(cos(yaw_old) - cos(yaw_old + yaw_rate*delta_t));;
			yaw_new = yaw_old + yaw_rate*delta_t;
		}
		else
		{
			x_new = x_old + velocity * delta_t * cos(yaw_old);
			y_new = y_old + velocity * delta_t * sin(yaw_old);
			yaw_new = yaw_old;
		}

		normal_distribution<double> dist_x(x_new, std_pos[0]);
		normal_distribution<double> dist_y(y_new, std_pos[1]);
		normal_distribution<double> dist_yaw(yaw_new, std_pos[2]);

		this->particles[i].x = dist_x(gen);
		this->particles[i].y = dist_y(gen);
		this->particles[i].theta = dist_yaw(gen);
		//agian, don't use this->particles[i].x += dist_x(gen), it confused me a long time before.
		this->particles[i].weight = 1.0;
		// every time we predict the next position, always set weight to 1.0
	}

}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	

	//observations vector is observation data you got from sensor ( such as laser )
	//predicted vector is landmark position you can find in the sensor range ( I still don't understand why called it as "predicted"... )
	std::vector<double> predicted_distance_vector;
	std::vector<int> predicted_id_vector;
	int minimal_distance_index;
	for (int i = 0; i < int(observations.size()); i++)
	{
		for (int j = 0; j < int(predicted.size()); j++)
		{
			predicted_distance_vector.push_back(dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y ));
			predicted_id_vector.push_back(predicted[j].id);
		}
		minimal_distance_index = std::distance(predicted_distance_vector.begin(), std::min_element(predicted_distance_vector.begin(), predicted_distance_vector.end())) ;
		observations[i].id = predicted_id_vector[minimal_distance_index];//for every observation we set a landmark id which is nearst to it, as it's id. 
		predicted_distance_vector.clear();
		predicted_id_vector.clear();
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

	
	std::vector<LandmarkObs> predicted;
	LandmarkObs pt;

	std::vector<LandmarkObs> transformed_observations;
	LandmarkObs ot;

	double x_part;
	double y_part;
	double theta;

	double x_obs;
	double y_obs;

	double x_map;
	double y_map;

	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];;

	double x_obs_transformed;
	double y_obs_transformed;

	double mu_x;
	double mu_y;

	double gauss_norm;
	double exponent;

	for (int i = 0; i < this->num_particles; i++)
	{
		//==================== part 1.for every particle, record those landmarks which are in senor range   ========================//
		predicted.clear();
		for (int j = 0; j <  int(map_landmarks.landmark_list.size()); j++)
		{
			if (dist(this->particles[i].x, this->particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range)
			{
				pt.id = map_landmarks.landmark_list[j].id_i;
				pt.x = map_landmarks.landmark_list[j].x_f;
				pt.y = map_landmarks.landmark_list[j].y_f;
				predicted.push_back(pt);
			}
		}
		//==================== part 1.for every particle, record those landmarks which are in senor range  ========================//

		//==================== part 2. coordinate transformation ========================//
		//transform observation from VEHICLE'S coordinate system to MAP'S coordinate system
		x_part = this->particles[i].x;
		y_part = this->particles[i].y;
		theta = this->particles[i].theta;

		transformed_observations.clear();
		//use a new particle vector to restore transformed observations
		for (int j = 0; j < int(observations.size()); j++)
		{
			x_obs = observations[j].x;
			y_obs = observations[j].y;
			x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
			y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

			ot.id = observations[j].id;
			ot.x = x_map;
			ot.y = y_map;
			transformed_observations.push_back(ot);
		}
		//==================== part 2. coordinate transformation ========================//

		//==================== part 3. dataAssociation ========================//
		dataAssociation(predicted , transformed_observations);
		//==================== part 3. dataAssociation ========================//
		
		//==================== part 4. calculate weight for every observation and multiply all weights ========================//
		this->particles[i].weight = 1.0;
		// make sure initial weight of every particle == 0;
		for (int j = 0; j < int(transformed_observations.size()); j++)
		{

			x_obs_transformed = transformed_observations[j].x;
			y_obs_transformed = transformed_observations[j].y;
			for (int n = 0; n < int(predicted.size()); n++)
			{
				if (predicted[n].id == transformed_observations[j].id)
				{
					mu_x = predicted[n].x;
					mu_y = predicted[n].y;
					break;
				}
			}
			gauss_norm = (1.0 / (2.0 * M_PI * sig_x * sig_y));
			exponent = ((x_obs_transformed - mu_x)*(x_obs_transformed - mu_x)) / (2 * sig_x * sig_x) + ((y_obs_transformed - mu_y) * (y_obs_transformed - mu_y)) / (2.0 * sig_y * sig_y);
			this->particles[i].weight = this->particles[i].weight * gauss_norm * exp(-exponent);
		}
		//==================== part 4. calculate weight for every observation and multiply all weights ========================//


		this->weights.push_back(this->particles[i].weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<Particle> resample_particles;
	std::vector<double> weights;
	for (int i = 0; i < this->num_particles; i++)
	{
		weights.push_back(this->particles[i].weight);
	}
	
	std::discrete_distribution<> dist(std::begin(weights), std::end(weights));
	std::random_device seed;
	std::mt19937 random_generator(seed());

	for (int i = 0; i < this->num_particles; i++)
	{
		resample_particles.push_back(this->particles[dist(random_generator)]);
	}
	this->particles = resample_particles;

	//reference:https://discussions.udacity.com/t/output-always-zero/260432/
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
