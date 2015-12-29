//==================================================================================================
// Name        : simpleMPM.h
// Author      : Ganzenmueller/Sandoval M.
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Basic Material Point Method Algorithm in C++, Ansi-style
// Created     : 15.12.2015
//==================================================================================================

#ifndef SIMPLEMPM_H_
#define SIMPLEMPM_H_

#include <Eigen/Dense>
using namespace Eigen;

struct Gridnode {

	double   m;                   // mass at grid node
	Vector2d q, v, a, fInt, fExt; // momentum, velocity, acceleration, forces.

	Gridnode() { // default initialization

		m = 0.0;
		q.setZero();
		v.setZero();
		a.setZero();
		fInt.setZero();
		fExt.setZero();
	}
};

struct Particle {
	double   m, V;             // particle mass
	Vector2d x, v, a;          // particle position and velocities
	Vector2d v_ipol, acc_ipol; // interpolated (from grid) new velocities and accelerations
	Matrix2d L;
	Matrix2d F;// velocity gradient tensor
	Matrix2d EPS;              // strain tensor
	Matrix2d SIG;              // stress tensor

	Particle() { // initialize particle data to some default values
		m = 0.0;
		V = 0.0;
		x.setZero();
		v.setZero();
		a.setZero();
		L.setZero();
		F.setZero();
		EPS.setZero();
		SIG.setZero();
	}
};

int npx;
int npy;
int Np;  // number of particles

int nnx; // number of grid nodes along x
int nny; // number of grid nodes along y
int Nn;  // total number of cells

int nTimeStep;
double p_spacing, cellSize, icellSize, dt; // discretization parameters
double xmin, xmax, ymin, ymax;
double rho, K, G;
Vector2d gravity;

Particle *particles;
Gridnode *gridnodes;

void   matProperties();
void   dirichletBC();
double weightFunction(double);
double derivativeWeightF(double);

#endif /* SIMPLEMPM_H_ */
