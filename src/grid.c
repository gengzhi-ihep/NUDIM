#include "grid.h"

/* rotation axes is along z, X is horizantal and Y is vertical.
 * Note x stores a slice with dimension nx*nz */
static void grid_slice_3D(int nx, int nz, double gama, double *x){

	double intervalx, intervalz, x0, z0;
	int i, j;

	intervalx = 1.0 / nx;
	intervalz = 1.0 / nz;

	z0 = -0.5-intervalz;

	for(j=0; j<nz; j++){

		z0 += intervalz;

		x0 = -0.5-intervalx;

		for(i=0; i<nx; i++){

			x0 += intervalx;

			x[3*(j*nx+i)+0] = x0 * cos(gama/180*PI);
			x[3*(j*nx+i)+1] = x0 * sin(gama/180*PI);
			x[3*(j*nx+i)+2] = z0;

		}
	}

	return;

}


extern double* equal_angle_sample_3D(double angle_min, double angle_max, double angle_interval,
									 int nx, int nz, int *ntotal){

	double *x, angle;
	int nslices, i, j;

	nslices = (int) ( (angle_max - angle_min)/angle_interval + 0.5 ) + 1;

	*ntotal = nslices * nx * nz;
	x = (double*) malloc(nslices*nx*nz*3*sizeof(double));

	angle = angle_min - angle_interval;

	for(i=0; i<nslices; i++){
		angle += angle_interval;
		grid_slice_3D(nx, nz, angle, &(x[i*nx*nz*3]));
	}

	return x;

}

/* rotation axes is along z, X is horizantal and Y is vertical
 * Note x stores a line with dimension nx. */
static void grid_slice_2D(int nx, double gama, double *x){

	double interval, x0;
	int i;

	interval = 1.0 / nx;
	x0 = -0.5 - interval;

	for(i=0; i<nx; i++){
		x0 += interval;
		x[2*i+0] = x0 * sin(gama/180*PI);
		x[2*i+1] = x0 * cos(gama/180*PI);
	}

	return;
}

extern double* equal_angle_sample_2D(double angle_min, double angle_max, double angle_interval,
									 int nx, int *ntotal){

	double *x, angle;
	int nslices, i;

	nslices = (int) ( (angle_max - angle_min)/angle_interval + 0.5 ) + 1;

	*ntotal = nslices * nx;
	x = (double*) malloc(nslices*nx*2*sizeof(double));

	angle = angle_min - angle_interval;

	for(i=0; i<nslices; i++){
		angle += angle_interval;
		grid_slice_2D(nx, angle, &(x[i*nx*2]));
	}

	return x;
}

/* uniform sampling 2D grid */
extern double* uniform_sample_2D(int nx, int nz, int *ntotal){

	double *x;
	double intervalx, intervalz, x0, z0;
	int i, j;

	x = (double*) malloc(nx*nz*2*sizeof(double));

	*ntotal = nx*nz;
	intervalx = 1.0 / nx;
	intervalz = 1.0 / nz;

	z0 = -0.5 -intervalz;

	for(j=0; j<nz; j++){

		z0 += intervalz;

		x0 = -0.5 - intervalx;

		for(i=0; i<nx; i++){

			x0 += intervalx;

			x[2*(j*nx+i)+0] = x0;
			x[2*(j*nx+i)+1] = z0;

		}

	}

	return x;

}

