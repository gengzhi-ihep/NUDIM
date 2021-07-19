#ifndef __DENOISE__

#define __DENOISE__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double ***total_variation(double*** img, float lam, float dt, int niter,
						 int nx, int ny, int nz);

double ***total_variation_discrete(double*** img, float lam, int niter,
						           int nx, int ny, int nz);

double ***total_variation_discrete_height(double*** img, float lam, int niter,
						           int nx, int ny, int nz, int nheight, int ncenter);

double ***aniso_diffusion(double*** img, float lambda,float kappa, int niter,
						 int nx, int ny, int nz);

double ***aniso_diffusion_2d(double*** img, float lambda,float kappa, int niter,
						 int nx, int ny, int nz);

double ***show_grad_mag(double***img, int nx, int ny, int nz);

double ***malloc_image(int nx, int ny, int nz);

void free_image(double*** img, int nx, int ny);

#endif
