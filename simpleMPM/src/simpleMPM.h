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
	Vector2d q, v, dv, a, fInt, fExt; // momentum, velocity, acceleration, forces.

	Gridnode() { // default initialization

		m = 0.0;
		q.setZero();
		v.setZero();
		dv.setZero();
		a.setZero();
		fInt.setZero();
		fExt.setZero();
	}
};

struct Particle {
	int      typ;
	double   m, V0, V;             // particle mass
	Vector2d agrid, dvgrid, vgrid, v, x ;          // particle position and velocities
	Vector2d v_ipol, acc_ipol; // interpolated (from grid) new velocities and accelerations
	Matrix2d L;
	Matrix2d F;// velocity gradient tensor
	Matrix2d EPS;              // strain tensor
	Matrix2d SIG;              // stress tensor

	Particle() { // initialize particle data to some default values
		typ = 0;
		m   = 0.0;
		V0  = 0.0;
		V   = 0.0;
		x.setZero();
		v.setZero();
		vgrid.setZero();
		dvgrid.setZero();
		agrid.setZero();
		L.setZero();
		F.setIdentity();
		EPS.setZero();
		SIG.setZero();
	}
};

struct IsoMat {

	double rho, K, G, nu, cp, k;

	IsoMat() {

		rho = 1.0;
		K   = 1.0;
		G   = 1.0;
		nu  = 1.0;
		cp  = 1.0;
		k   = 1.0;

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
Vector2d gravity;

Particle *particles;
Gridnode *gridnodes;
IsoMat   *materials;

void   matProperties();
void   dirichletBC();
double weightFunction(double);
double derivativeWeightF(double);
Matrix2d Deviator(const Matrix2d M);

#endif /* SIMPLEMPM_H_ */
