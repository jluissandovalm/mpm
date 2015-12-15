//============================================================================
// Name        : simpleMPM.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "simpleMPM.h"
#include <iostream>
#include <stdio.h>
using namespace std;

#define BIG 1.0e22;

void settings() {
	/*
	 * define global simulation box settings
	 */

	Np = 2; // number of particles
	particles = new Particle[Np];
	cellsize = 1.0;
	icellsize = 1.0 / cellsize;
}

void createParticles() {
	/*
	 * create initial particle layout
	 */

	// two particles
	particles[0].x(0) = 10.1;
	particles[0].x(1) = -1.0e-16;
	particles[0].m = 1.0;
	particles[1].x(0) = 20.1;
	particles[1].m = 1.0;


}

void createGrid() {

	// 1) find minimum and maximum particle coordinates
	xmin = ymin = BIG
	;
	xmax = ymax = -BIG
	;
	for (int i = 0; i < Np; i++) {
		if (particles[i].x(0) > xmax)
			xmax = particles[i].x(0);
		if (particles[i].x(0) < xmin)
			xmin = particles[i].x(0);
		if (particles[i].x(1) > ymax)
			ymax = particles[i].x(1);
		if (particles[i].x(1) < ymin)
			ymin = particles[i].x(1);
	}

	xmin = floor(xmin); // this ensures the grid does not move with the particles
	ymin = floor(ymin);

	xmin -= cellsize; // this ensures we have only positive cell indices
	ymin -= cellsize;
	xmax -= cellsize;
	ymax -= cellsize;

	printf("particles coordinates: %f < x < %f\n", xmin, xmax);
	printf("particles coordinates: %f < y < %f\n", ymin, ymax);

	// 2) find number of grid nodes along x and y
	Nx = icellsize * (xmax - xmin) + 4; // need two extra cells on either side
	Ny = icellsize * (ymax - ymin) + 4;
	Ncells = Nx * Ny;

	// allocate storage for grid
	gridnodes = new Gridnode[Ncells];

}

void pointsToGrid() {
	/*
	 * This routine implements Buzzi's "Caveat ..." method for iterating over all grid nodes which are visible from a given particle.
	 *
	 * Here, you need to interpolate particle momentum and mass to grid nodes,
	 */
	int ref_node, node_index, ix, iy;
	int ref_ix, ref_iy;
	double dx, dy, wfx, wfy, wf;

	for (int i = 0; i < Np; i++) {

		ref_ix = (int) (icellsize * (particles[i].x(0) - xmin)) - 1; // this is the x location of the reference node in cell space
		ref_iy = (int) (icellsize * (particles[i].x(1) - ymin)) - 1; // this is the y location of the reference node in cell space
		ref_node = ref_ix + ref_iy * Nx; // this is the index of the reference node

		for (ix = 0; ix < 4; ix++) { // this is the x stencil visible to the current particle

			// x distance:
			dx = cellsize * (ref_ix + ix) - (particles[i].x(0) - xmin); // this is the x distance in world coords
			printf("particle x =%f, grid x =%f, distance =%f\n", particles[i].x(0) - xmin, cellsize * (ref_ix + ix), dx);
			// compute weight function for dx
			wfx = 0.0;

			for (iy = 0; iy < 4; iy++) { // this is the y stencil visible to the current particle

				// y distance:
				dy = cellsize * (ref_iy + iy) - (particles[i].x(1) - ymin);
				printf("particle y =%f, grid y =%f, distance =%f\n", particles[i].x(1) - ymin, cellsize * (ref_iy + iy), dy);
				wfy = 0.0; // compute weight function for dy

				node_index = ref_node + ix + iy * Nx; // current node index

				if ((node_index < 0) || (node_index > Ncells - 1)) { // memory access error check
					printf("node index %d outside allowed range %d to %d\n", node_index, 0, Ncells);
					exit(1);
				}

				//printf("particle %d has ref node %d, node_index %d\n", i, ref_node, node_index);
				wf = wfx * wfy; // total weight function dyadic product;

				// interpolate particle values to node
				gridnodes[node_index].m += wf * particles[i].m;

			} // end loop over y dimension on grid
		} // end loop over x dimension on grid
	} // end loop over particles
}

void computeParticleGradients() {
	/*
	 * iterate over all particles and compute the grid field gradient at the particle locations,
	 * e.g.: compute the gradient of the velocity field
	 */

	// loop over all particles
	// loop over all grid points visible to a particle
	// compute weights and weight derivatives
	// compute velocity gradient
}

void computeGridRates() {
	/*
	 * determine time derivative of grid field.
	 * for now
	 * 1) compute a particle stress rate using the velocity gradient,
	 * 2) update the particle stress state using its stress rate
	 * 3) iterate over all visible grid nodes and compute a force on the grid nodes from the divergence of the stress field
	 */
}

void advanceGrid() {
	/*
	 * given the grid forces, update the grid velocities
	 */
}

void gridToParticles() {
	/*
	 * loop over all particles and their respective visible grid nodes and compute
	 * 1) interpolated grid velocities at particle location
	 * 2) interpolated grid accelerations at particle location
	 */
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
	 * 	1) with the interpolated grid velocities at the particle location, update particle positions
	 * 	2) with the interpolated grid accelerations at the particle location, update particle velocity.
	 */
}

void dumpParticles() {
	/*
	 * dump particle configuration to file format readable by OVITO
	 */
}

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	settings();
	createParticles();

	// loop over timesteps
	for (int ntimestep = 0; ntimestep < 1; ntimestep++) {
		createGrid();
		pointsToGrid();
		computeParticleGradients();
		computeGridRates(); // grid node accelerations
		advanceGrid(); // time integration of grid
		gridToParticles(); // gather values from grid
		destroyGrid();

		advanceParticles();
		dumpParticles();
	}

	return 0;
}
