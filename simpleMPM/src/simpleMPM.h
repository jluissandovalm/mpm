/*
 * simpleMPM.h
 *
 *  Created on: Dec 15, 2015
 *      Author: ganzenmueller
 */

#ifndef SIMPLEMPM_H_
#define SIMPLEMPM_H_

#include <Eigen/Dense>
using namespace Eigen;

struct Gridnode {
	Vector2d v, f; // grid velocities and forces
	double m; // mass at grid node

	Gridnode() { // default initialization
		v.setZero();
		f.setZero();
		m = 0.0;
	}

};

struct Particle {
	Vector2d x, v; // particle position and velocities
	Vector2d v_ipol, acc_ipol; // interpolated (from grid) new velocities and accelerations
	double m; // particle mass

	Particle() { // initialize particle data to some default values
		x.setZero();
		v.setZero();
		m = 0.0;
	}
};

int Ncells; // total number of cells
int Np; // number of particles
int Nx; // number of grid nodes along x
int Ny; // number of grid nodes along y
double cellsize, icellsize; // discretization width and its inverse
double xmin, xmax, ymin, ymax;
Particle *particles;
Gridnode *gridnodes;



#endif /* SIMPLEMPM_H_ */
