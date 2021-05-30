#ifndef __GRID__

#define __GRID__

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926

static void grid_slice_3D(int nx, int nz, double gama, double *x);

extern double* equal_angle_sample_3D(double angle_min, double angle_max, double angle_interval,
									 int nx, int nz, int *ntotal);

static void grid_slice_2D(int nx, double gama, double *x);

extern double* equal_angle_sample_2D(double angle_min, double angle_max, double angle_interval,
									 int nx, int *ntotal);

extern double* uniform_sample_2D(int nx, int nz, int *ntotal);

#endif
