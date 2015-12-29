//==================================================================================================
// Name        : simpleMPM.cpp
// Author      : Ganzenmüller/Sandoval M.
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : Basic Material Point Method Algorithm in C++, Ansi-style
// Created     : 15.12.2015
//==================================================================================================

#include "simpleMPM.h"
#include <iostream>
#include <stdio.h>
using namespace std;

#define BIG 1.0e22
#define MASS_CUTOFF 1.0E-16

/********************************** M A I N   S U B R U T I N E S *********************************/

void settings() {
	/*
	 * define global simulation box settings
	 */

	double c0, sf;

	npx = 20;
	npy = 10;
	Np  = npx * npy;

	p_spacing = 1.0;
	cellSize  = 1.0 * p_spacing;
	icellSize = 1.0 / cellSize;

	matProperties();

	c0 = sqrt(K / rho);
	sf = 0.80;

	dt = sf * p_spacing / c0;
	printf("initial_dt = %f\n", dt);

	nTimeStep = 0;

	gravity << 0.000, -0.000;

}

void createParticles() {
	/*
	 * create initial particle layout
	 */

	particles = new Particle[Np];

	int p_index = 0;

	for (int j = 0; j < npy; j++) {
		for (int i = 0; i < 10; i++) {

			particles[p_index].x(0) = i * p_spacing + 0.5;
			particles[p_index].x(1) = j * p_spacing + 1.5;
			particles[p_index].v(0) = +0.80;
			particles[p_index].v(1) =  0.00;
			particles[p_index].m    = 1.0;
			particles[p_index].V    = p_spacing * p_spacing * p_spacing;

			p_index++;
		}
	}

	for (int j =  0; j < npy; j++) {
		for (int i = 10; i < npx; i++) {
			particles[p_index].x(0) = i * p_spacing + 4.5;
			particles[p_index].x(1) = j * p_spacing + 1.5;
			particles[p_index].v(0) = -0.80;
			particles[p_index].v(1) =  0.00;
			particles[p_index].m    = 1.0;
			particles[p_index].V    = p_spacing * p_spacing * p_spacing;

			p_index++;
		}
	}

}


void createGrid() {

	// 1) find minimum and maximum particle coordinates

	xmin = ymin =  BIG;
	xmax = ymax = -BIG;

	for (int p_index = 0; p_index < Np; p_index++) {

		if (particles[p_index].x(0) > xmax)   xmax = particles[p_index].x(0);

		if (particles[p_index].x(0) < xmin)   xmin = particles[p_index].x(0);

		if (particles[p_index].x(1) > ymax)   ymax = particles[p_index].x(1);

		if (particles[p_index].x(1) < ymin)   ymin = particles[p_index].x(1);
	}

	// this ensures the grid does not move with the particles
	xmin = floor(xmin);
	ymin = floor(ymin);
	xmax = floor(xmax); // ?
	ymax = floor(ymax); // ?

	// this ensures we have only positive cell indices
	xmin -= cellSize;
	ymin -= cellSize;
	xmax -= cellSize;
	ymax -= cellSize;

	printf("particles coordinates:\n%f < x < %f\n%f < y < %f\n", xmin, xmax, ymin, ymax);

	// 2) find number of grid nodes along x and y
	nnx = icellSize * (xmax - xmin) + 4; // "... + 4" because need two extra cells on either side
	nny = icellSize * (ymax - ymin) + 4;
	Nn  = nnx * nny;

	// allocate storage for grid
	gridnodes = new Gridnode[Nn];

}


void particlesToGrid() {
	/*
	 * This routine implements Buzzi's "Caveat ..." method for iterating over all grid nodes which
	 * are visible from a given particle.
	 *
	 * Here, you need to interpolate mass and particle momentum to grid nodes,
	 */
	int n_ref, n_index, ix, iy;
	int ref_ix, ref_iy;
	double dx, dy, wfx, wfy, wf;

	for (int p_index = 0; p_index < Np; p_index++) {       // Loop over all particles

		// this is the x location of the reference node in cell space
		ref_ix = (int) (icellSize * (particles[p_index].x(0) - xmin)) - 1;
		// this is the y location of the reference node in cell space
		ref_iy = (int) (icellSize * (particles[p_index].x(1) - ymin)) - 1;
		// this is the index of the reference node
		n_ref = ref_ix + ref_iy * nnx;

		for (ix = 0; ix < 4; ix++) {     // Loop over x-stencil visible to the current particle

			// x distance in grid coords:
			dx = (ref_ix + ix) - icellSize * (particles[p_index].x(0) - xmin);
			// printf("particle x = %f, grid x = %f, distance = %f\n",
			//        particles[i].x(0) - xmin, cellSize * (ref_ix + ix), dx);

			// compute weight function for dx
			wfx = weightFunction(dx);

			for (iy = 0; iy < 4; iy++) { // Loop over y-stencil visible to the current particle

				// y distance in grid coords:
				dy = (ref_iy + iy) - icellSize * (particles[p_index].x(1) - ymin);
				// printf("particle y = %f, grid y = %f, distance = %f\n",
				//        particles[i].x(1) - ymin, cellSize * (ref_iy + iy), dy);

				// compute weight function for dy
				wfy = weightFunction(dy);

				n_index = n_ref + ix + iy * nnx; // current node index

				if ((n_index < 0) || (n_index > Nn - 1)) { // memory access error check
					printf("node index %d outside allowed range %d to %d\n", n_index, 0, Nn);
					exit(1);
				}

				wf = wfx * wfy; // total weight function (kernel) dyadic product

				// interpolate particle values to node
				gridnodes[n_index].m += wf * particles[p_index].m;
				gridnodes[n_index].q += wf * particles[p_index].m * particles[p_index].v;

			} // end loop over y dimension on grid
		} // end loop over x dimension on grid
	} // end loop over particles

	//dirichletBC();

	for (int n_index = 0; n_index < Nn; n_index++) {

		if (gridnodes[n_index].m > MASS_CUTOFF) {

			gridnodes[n_index].v = gridnodes[n_index].q / gridnodes[n_index].m;

		}
	}

}


void computeParticleGradients() {

	/* Iterate over all particles and compute the grid field gradient at the particle locations,
	 * e.g.: compute the gradient of the velocity field
	 *
	 * 1) loop over all particles
	 * 2) loop over all grid points visible to a particle
	 * 3) compute weights and weight derivatives
	 * 4) compute velocity gradient
	 */

	int n_ref, n_index, ix, iy;
	int ref_ix, ref_iy;
	double dx, dy, wfx, wfy, dwfx, dwfy;
	Vector2d gwf;

	for (int p_index = 0; p_index < Np; p_index++) {

		ref_ix = (int) (icellSize * (particles[p_index].x(0) - xmin)) - 1;
		ref_iy = (int) (icellSize * (particles[p_index].x(1) - ymin)) - 1;
		n_ref = ref_ix + ref_iy * nnx;

		for (ix = 0; ix < 4; ix++) {

			dx = (ref_ix + ix) - icellSize * (particles[p_index].x(0) - xmin);

			wfx = weightFunction(dx);
			dwfx = derivativeWeightF(dx);

			for (iy = 0; iy < 4; iy++) {

				dy = (ref_iy + iy) - icellSize * (particles[p_index].x(1) - ymin);

				wfy = weightFunction(dy);
				dwfy = derivativeWeightF(dy);

				n_index = n_ref + ix + iy * nnx;

				if ((n_index < 0) || (n_index > Nn - 1)) {
					printf("node index %d outside allowed range %d to %d\n", n_index, 0, Nn);
					exit(1);
				}

				gwf << dwfx * wfy, wfx * dwfy;
				particles[p_index].L   += gridnodes[n_index].v * gwf.transpose();

			} // end loop over y dimension on grid
		} // end loop over x dimension on grid
	} // end loop over particles

}


void computeGridRates() {

	/* Determine time derivative of grid field.
	 * for now
	 * 1) Compute a particle stress rate using the velocity gradient.
	 * 2) Update the particle stress state using its stress rate.
	 * 3) Iterate over all visible grid nodes and compute a force on the grid nodes from the
	 *    divergence of the stress field.
	 */

	double Vnew, pressure;
	Matrix2d D, W, devD, devE, II;
	II.setIdentity();

	for (int p_index = 0; p_index < Np; p_index++) {

		D = 0.5 * (particles[p_index].L + particles[p_index].L.transpose()); // stretch rate tensor
		W = 0.5 * (particles[p_index].L - particles[p_index].L.transpose()); // spin    rate tensor


		particles[p_index].EPS += dt * D;
		devE = particles[p_index].EPS - (1.0/2.0) * particles[p_index].EPS.trace() * II;
		particles[p_index].SIG  = K * particles[p_index].EPS.trace() * II + 2 * G * devE;

		/*
		particles[p_index].F += dt * particles[p_index].L * particles[p_index].F;
		Vnew = particles[p_index].V * particles[p_index].F.determinant();
		*/
		/*
		devD = 2.0 * G * (D - (1.0/2.0) * D.trace() * II);
		Vnew = particles[p_index].V + dt * particles[p_index].V * D.trace();
		pressure = K * (particles[p_index].V / Vnew - 1.0);
		particles[p_index].V   =  Vnew;
		particles[p_index].SIG = devD - pressure * II;
		*/
	}


	int n_ref, n_index, ix, iy;
	int ref_ix, ref_iy;
	double dx, dy, wfx, wfy, dwfx, dwfy, wf;
	Vector2d gwf;

	for (int p_index = 0; p_index < Np; p_index++) {

			ref_ix = (int) (icellSize * (particles[p_index].x(0) - xmin)) - 1;
			ref_iy = (int) (icellSize * (particles[p_index].x(1) - ymin)) - 1;
			n_ref = ref_ix + ref_iy * nnx;

			for (ix = 0; ix < 4; ix++) {

				dx = (ref_ix + ix) - icellSize * (particles[p_index].x(0) - xmin);

				wfx = weightFunction(dx);
				dwfx = derivativeWeightF(dx);

				for (iy = 0; iy < 4; iy++) {

					dy = (ref_iy + iy) - icellSize * (particles[p_index].x(1) - ymin);

					wfy = weightFunction(dy);
					dwfy = derivativeWeightF(dy);

					n_index = n_ref + ix + iy * nnx;

					if ( (n_index < 0) || (n_index > Nn - 1) ) {
						printf("node index %d outside allowed range %d to %d\n", n_index, 0, Nn);
						exit(1);
					}

					wf   = wfx * wfy;
					gwf << dwfx * wfy, wfx * dwfy;

					gridnodes[n_index].fInt += particles[p_index].V * particles[p_index].SIG * gwf;
					gridnodes[n_index].fExt += particles[p_index].m * gravity * wf /* + NBC*/;

				} // end loop over y dimension on grid
			} // end loop over x dimension on grid
		} // end loop over particles

}


void advanceGrid() {

	/* Given the grid forces, update the grid velocities
	 */

	Vector2d q_dot;

	for (int n_index = 0; n_index < Nn; n_index++) {

		q_dot = gridnodes[n_index].fExt - gridnodes[n_index].fInt;

		gridnodes[n_index].q +=                    q_dot * dt;
		gridnodes[n_index].a  =                    q_dot / gridnodes[n_index].m;
		gridnodes[n_index].v  = gridnodes[n_index].q     / gridnodes[n_index].m;

	}

	//dirichletBC();

}

void gridToParticles() {

	/* Loop over all particles and their respective visible grid nodes and compute
	 *
	 * 1) Interpolated grid velocities    at particle location
	 * 2) Interpolated grid accelerations at particle location
	 */

	int n_ref, n_index, ix, iy;
		int ref_ix, ref_iy;
		double dx, dy, wfx, wfy, wf;
		Vector2d gwf;

		for (int p_index = 0; p_index < Np; p_index++) {

				ref_ix = (int) (icellSize * (particles[p_index].x(0) - xmin)) - 1;
				ref_iy = (int) (icellSize * (particles[p_index].x(1) - ymin)) - 1;
				n_ref = ref_ix + ref_iy * nnx;

				for (ix = 0; ix < 4; ix++) {

					dx = (ref_ix + ix) - icellSize * (particles[p_index].x(0) - xmin);

					wfx = weightFunction(dx);

					for (iy = 0; iy < 4; iy++) {

						dy = (ref_iy + iy) - icellSize * (particles[p_index].x(1) - ymin);

						wfy = weightFunction(dy);

						n_index = n_ref + ix + iy * nnx;

						if ( (n_index < 0) || (n_index > Nn - 1) ) {
							printf("node index %d outside allowed range %d to %d\n", n_index, 0, Nn);
							exit(1);
						}

						wf   = wfx * wfy;

						particles[p_index].v += dt * wf * gridnodes[n_index].a;
						particles[p_index].a +=      wf * gridnodes[n_index].a;

					} // end loop over y dimension on grid
				} // end loop over x dimension on grid
			} // end loop over particles

}


void destroyGrid() {
	/*
	 * we delete the grid at the end of this time step
	 */
	delete[] gridnodes;
}


void advanceParticles() {
	/*
	 * loop over all particles
	 * 	1) with the interpolated grid velocities,    update particle positions
	 * 	2) with the interpolated grid accelerations, update particle velocity.
	 */

	for (int p_index = 0; p_index < Np; p_index++) {

		particles[p_index].x += dt * particles[p_index].v;

	}

}


void dumpGrid() {
	/*
	 * dump particle configuration to file format readable by OVITO
	 */

	FILE *fp = fopen("grid.dump", "a");
	fprintf(fp, "ITEM: TIMESTEP\n%d\n", nTimeStep);
	fprintf(fp, "ITEM: NUMBER OF ATOMS\n%d\n", Nn + Np);
	fprintf(fp, "ITEM: BOX BOUNDS\n  %8.4f %8.4f\n  %8.4f %8.4f\n  %8.4f %8.4f\n",
			xmin, xmax + 3 * cellSize, ymin, ymax + 3 * cellSize, -0.1, 0.1);
	fprintf(fp,	"ITEM: ATOMS ID TYPE X Y Z Vx Vy Vz Mass Vol Lxx Lxy Lyx Lyy SIGxx SIGxy SIGyx SIGyy\n");

	for (int iy = 0; iy < nny; iy++) {
		for (int ix = 0; ix < nnx; ix++) {
			int n_index = ix + iy * nnx;

			if ( (n_index < 0) || (n_index > Nn - 1) ) {
				printf("node index %d outside allowed range %d to %d\n", n_index, 0, Nn);
				exit(1);
			}

			double xcoord = ix * cellSize + xmin;
			double ycoord = iy * cellSize + ymin;

			fprintf(fp, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", n_index + 1, 1,
					xcoord,                   ycoord,                  0.0,
					gridnodes[n_index].v(0),  gridnodes[n_index].v(1), 0.0,
					gridnodes[n_index].m,     0.0,
					0.0, 0.0, 0.0, 0.0,
					0.0, 0.0, 0.0, 0.0);
		}

	}

	for (int i = 0; i < Np; i++) {

		fprintf(fp, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", i + 10000, 2,
				particles[i].x(0),     particles[i].x(1),     1.0 * p_spacing,
		        particles[i].v(0),     particles[i].v(1),     0.0,
				particles[i].m,        particles[i].V,
				particles[i].L(0,0),   particles[i].L(0,1),
				particles[i].L(1,0),   particles[i].L(1,1),
				particles[i].SIG(0,0), particles[i].SIG(0,1),
				particles[i].SIG(1,0), particles[i].SIG(1,1)
				);
	}
	fclose(fp);

}


void dumpParticles() {
	/*
	 * dump particle configuration to file format readable by OVITO
	 */

}

/****************************** A U X I L I A R Y   F U N C T I O N S *****************************/

void matProperties() {

	rho = 1.35E-5;
	K   = 2.60;
	G   = 1.05;

}

double weightFunction(double dx) {

	double r = fabs(dx);

	if (r < 1.0) {
		return 0.5 * r * r * r - r * r + 2.0 / 3.0;
	} else if (r < 2.0) {
		return -(r * r * r / 6.0) + r * r - 2.0 * r + 4.0 / 3.0;
	} else {
		return 0.0;
	}
}


double derivativeWeightF(double dx) {

	double r = fabs(dx);
	double sign;
	if (dx < 0.) {
		sign = 1.0 * icellSize;
	} else {
		sign = -1.0 * icellSize;
	}

	if (r < 1.0) {
		return sign * (1.5 * r * r - 2.0 * r);
	} else if (r < 2.0) {
		return sign * (-(r * r / 2.0) + 2 * r - 2.0);
	} else {
		return 0.0;
	}
}


void dirichletBC() {

	for (int n_index = 0; n_index < 12; n_index++) {

		gridnodes[n_index + 2 + nnx].q(0) = 0;
		gridnodes[n_index + 2 + nnx].q(1) = 0;
		gridnodes[n_index + 2 + nnx].a(0) = 0;
		gridnodes[n_index + 2 + nnx].a(1) = 0;

	}

}


/******************************************* M  A  I  N *******************************************/

int main() {

	cout << "!!!Simple GMPM!!!" << endl;

	settings();
	matProperties();
	createParticles();
	remove("grid.dump");
	int dump;

	// loop over timesteps
	for (nTimeStep = 0; nTimeStep < 8000; nTimeStep++) {

		printf("Time Step: %d\n", nTimeStep);

		createGrid();
		particlesToGrid();
		computeParticleGradients();
		computeGridRates(); // grid node accelerations

		if ( dump == nTimeStep / 10) {
			printf("dump = %d\n", dump);
			dumpGrid();
			dump++;
		}

		advanceGrid();      // time integration of grid
		gridToParticles();  // gather values from grid
		destroyGrid();

		advanceParticles();
		dumpParticles();
	}

	return 0;
}
