#include "denoise.h"

double*** total_variation(double*** img, float lam, float dt, int niter,
					 int nx, int ny, int nz){

	double ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uxz, ux2, uy2, uz2;
	double ***denoised;
	double div_grad;
	double ep = 1, ep2;
	int kp, km, jp, jm, ip, im;
	int i, j, k;

	ep2 = pow(ep,2);

	denoised = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(denoised[i][j],img[i][j],sizeof(double)*nz);
	}

	while(niter--){

	for(k=0; k<nz; k++){

		kp = (k==nz-1) ? nz-1 : k+1;
		km = (k == 0) ? 0 : k-1;

		for(j=0; j<ny; j++){

			jp = (j==ny-1)? ny-1: j+1;
			jm = (j==0) ? 0 : j-1;

			for(i=0; i<nx; i++){

				ip = (i==nx-1)? nx-1: i+1;
				im = (i==0) ? 0 : i-1;

				ux = (denoised[ip][j][k]-denoised[im][j][k])/2;
				uy = (denoised[i][jp][k]-denoised[i][jm][k])/2;
				uz = (denoised[i][j][kp]-denoised[i][j][km])/2;
				uxx = denoised[ip][j][k]+denoised[im][j][k]-2*denoised[i][j][k];
				uyy = denoised[i][jp][k]+denoised[i][jm][k]-2*denoised[i][j][k];
				uzz = denoised[i][j][kp]+denoised[i][j][km]-2*denoised[i][j][k];
				uxy = (denoised[ip][jp][k]+denoised[im][jm][k]-denoised[im][jp][k]-denoised[ip][jm][k])/4;
				uyz = (denoised[i][jp][kp]+denoised[i][jm][km]-denoised[i][jp][km]-denoised[i][jm][kp])/4;
				uxz = (denoised[ip][j][kp]+denoised[im][j][km]-denoised[ip][j][km]-denoised[im][j][kp])/4;
				ux2 = pow(ux,2.0);
				uy2 = pow(uy,2.0);
				uz2 = pow(uz,2.0);
				div_grad = ((ux2+ep2)*(uyy+uzz)+(uxx+uzz)*(uy2+ep2)+(uz2+ep2)*(uxx+uyy)
							-2*(ux*uy*uxy+ux*uz*uxz+uy*uz*uyz))/pow((ux2+uy2+uz2+ep2),1.5);

				denoised[i][j][k] += dt*(div_grad - lam*(denoised[i][j][k]-img[i][j][k]));

			}

		}
	}

	}

	return denoised;

}

double*** total_variation_discrete(double*** img, float lam, int niter, int nx, int ny, int nz){

	double ***denoised, ***grad;
	double eps = 1e-8, dA;
	double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12;
	int kp, km, jp, jm, ip, im;
	int i, j, k;

	denoised = malloc_image(nx, ny, nz);
	grad = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(denoised[i][j],img[i][j],sizeof(double)*nz);
	}

	while(niter--){

	tmp0 = 0;

	for(k=0; k<nz; k++){

		kp = (k==nz-1) ? nz-1 : k+1;
		km = (k == 0) ? 0 : k-1;

		for(j=0; j<ny; j++){

			jp = (j==ny-1)? ny-1: j+1;
			jm = (j==0) ? 0 : j-1;

			for(i=0; i<nx; i++){

				ip = (i==nx-1)? nx-1: i+1;
				im = (i==0) ? 0 : i-1;

				tmp1 = denoised[i][j][k] - denoised[im][j][k];
				tmp2 = denoised[i][j][k] - denoised[i][jm][k];
				tmp3 = denoised[i][j][k] - denoised[i][j][km];
				tmp4 = denoised[ip][j][k] - denoised[i][j][k];
				tmp5 = denoised[ip][j][k] - denoised[ip][jm][k];
				tmp6 = denoised[ip][j][k] - denoised[ip][j][km];
				tmp7 = denoised[i][jp][k] - denoised[im][jp][k];
				tmp8 = denoised[i][jp][k] - denoised[i][j][k];
				tmp9 = denoised[i][jp][k] - denoised[i][jp][km];
				tmp10 = denoised[i][j][kp] - denoised[im][j][kp];
				tmp11 = denoised[i][j][kp] - denoised[i][jm][kp];
				tmp12 = denoised[i][j][kp] - denoised[i][j][k];

				//gradient with respect to each pixel
				grad[i][j][k] = (tmp1+tmp2+tmp3)/sqrt(eps+pow(tmp1,2)+pow(tmp2,2)+pow(tmp3,2))-
					            tmp4/sqrt(eps+pow(tmp4,2)+pow(tmp5,2)+pow(tmp6,2))-
								tmp8/sqrt(eps+pow(tmp7,2)+pow(tmp8,2)+pow(tmp9,2))-
								tmp12/sqrt(eps+pow(tmp10,2)+pow(tmp11,2)+pow(tmp12,2));
				tmp0 += pow(grad[i][j][k],2);
			}
		}
	}

	//scale gradient for each pixel
	dA = 0;
        tmp1 = sqrt(tmp0);
	for(k=0; k<nz; k++)
	for(j=0; j<ny; j++)
	for(i=0; i<nx; i++){
		grad[i][j][k] /= tmp1;
		if(denoised[i][j][k] < 0) dA += pow(denoised[i][j][k],2);
	}

	//update image using gradient descent methods
	tmp1 = sqrt(dA);
	for(k=0; k<nz; k++)
	for(j=0; j<ny; j++)
	for(i=0; i<nx; i++){
		denoised[i][j][k] -= lam*tmp1*grad[i][j][k];
//		denoised[i][j][k] -= lam*grad[i][j][k];
	}


	}

	free_image(grad, nx, ny);

	return denoised;

}

double*** total_variation_discrete_height(double*** img, float lam, int niter, int nx, int ny, int nz, int nheight, int ncenter){

	double ***denoised, ***grad;
	double eps = 1e-8, dA;
	double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12;
	int kp, km, jp, jm, ip, im;
	int i, j, k;

	denoised = malloc_image(nx, ny, nz);
	grad = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(denoised[i][j],img[i][j],sizeof(double)*nz);
	}

	while(niter--){

	tmp0 = 0;

	for(k=0; k<nz; k++){

		kp = (k==nz-1) ? nz-1 : k+1;
		km = (k == 0) ? 0 : k-1;

		for(j=0; j<ny; j++){

                        if(j<(2*ncenter-nheight)/2 || j>(2*ncenter+nheight)/2)continue;

			jp = (j==ny-1)? ny-1: j+1;
			jm = (j==0) ? 0 : j-1;

			for(i=0; i<nx; i++){

				ip = (i==nx-1)? nx-1: i+1;
				im = (i==0) ? 0 : i-1;

				tmp1 = denoised[i][j][k] - denoised[im][j][k];
				tmp2 = denoised[i][j][k] - denoised[i][jm][k];
				tmp3 = denoised[i][j][k] - denoised[i][j][km];
				tmp4 = denoised[ip][j][k] - denoised[i][j][k];
				tmp5 = denoised[ip][j][k] - denoised[ip][jm][k];
				tmp6 = denoised[ip][j][k] - denoised[ip][j][km];
				tmp7 = denoised[i][jp][k] - denoised[im][jp][k];
				tmp8 = denoised[i][jp][k] - denoised[i][j][k];
				tmp9 = denoised[i][jp][k] - denoised[i][jp][km];
				tmp10 = denoised[i][j][kp] - denoised[im][j][kp];
				tmp11 = denoised[i][j][kp] - denoised[i][jm][kp];
				tmp12 = denoised[i][j][kp] - denoised[i][j][k];

				//gradient with respect to each pixel
				grad[i][j][k] = (tmp1+tmp2+tmp3)/sqrt(eps+pow(tmp1,2)+pow(tmp2,2)+pow(tmp3,2))-
					            tmp4/sqrt(eps+pow(tmp4,2)+pow(tmp5,2)+pow(tmp6,2))-
								tmp8/sqrt(eps+pow(tmp7,2)+pow(tmp8,2)+pow(tmp9,2))-
								tmp12/sqrt(eps+pow(tmp10,2)+pow(tmp11,2)+pow(tmp12,2));
				tmp0 += pow(grad[i][j][k],2);
			}
		}
	}

	//scale gradient for each pixel
	dA = 0;
        tmp1 = sqrt(tmp0);
	for(k=0; k<nz; k++){
	    for(j=0; j<ny; j++){

                if(j<(2*ncenter-nheight)/2 || j>(2*ncenter+nheight)/2)continue;

	        for(i=0; i<nx; i++){
		    grad[i][j][k] /= tmp1;
		    if(denoised[i][j][k] < 0) dA += pow(denoised[i][j][k],2);
	        }
            }
        }

	//update image using gradient descent methods
	tmp1 = sqrt(dA);
	for(k=0; k<nz; k++){
	    for(j=0; j<ny; j++){

                if(j<(2*ncenter-nheight)/2 || j>(2*ncenter+nheight)/2)continue;

	        for(i=0; i<nx; i++){
		    denoised[i][j][k] -= lam*tmp1*grad[i][j][k];
//		    denoised[i][j][k] -= lam*grad[i][j][k];
	        }
            }
        }


	}

	free_image(grad, nx, ny);

	return denoised;

}

double*** aniso_diffusion(double*** img, float lambda, float kappa, int niter,
						 int nx, int ny, int nz){

	double ***denoised;
	double grad_N, grad_S, grad_E, grad_W, grad_F, grad_B;
	double cN, cS, cE, cW, cF, cB, Ksq;
	int kp, km, jp, jm, ip, im;
	int i, j, k;

	denoised = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(denoised[i][j],img[i][j],sizeof(double)*nz);
	}

	Ksq = pow(kappa, 2);

	while(niter--){

		for(k=0; k<nz; k++){

			kp = (k==nz-1)? nz-1 : k+1;
			km = (k==0) ? 0 : k-1;

			for(j=0; j<ny; j++){

				jp = (j==ny-1)? ny-1: j+1;
				jm = (j==0) ? 0 : j-1;

				for(i=0; i<nx; i++){

					ip = (i==nx-1)? nx-1: i+1;
					im = (i==0) ? 0 : i-1;

					grad_N = denoised[i][jm][k] - denoised[i][j][k];
					grad_S = denoised[i][jp][k] - denoised[i][j][k];
					grad_E = denoised[ip][j][k] - denoised[i][j][k];
					grad_W = denoised[im][j][k] - denoised[i][j][k];
					grad_F = denoised[i][j][kp] - denoised[i][j][k];
					grad_B = denoised[i][j][km] - denoised[i][j][k];

					cN = 1.0 / (1.0 + pow(grad_N, 2) / Ksq);
					cS = 1.0 / (1.0 + pow(grad_S, 2) / Ksq);
					cE = 1.0 / (1.0 + pow(grad_E, 2) / Ksq);
					cW = 1.0 / (1.0 + pow(grad_W, 2) / Ksq);
					cF = 1.0 / (1.0 + pow(grad_F, 2) / Ksq);
					cB = 1.0 / (1.0 + pow(grad_B, 2) / Ksq);
/*
					cN = exp( -pow(grad_N, 2) / Ksq );
					cS = exp( -pow(grad_S, 2) / Ksq );
					cE = exp( -pow(grad_E, 2) / Ksq );
					cW = exp( -pow(grad_W, 2) / Ksq );
					cF = exp( -pow(grad_F, 2) / Ksq );
					cB = exp( -pow(grad_B, 2) / Ksq );
*/
					denoised[i][j][k] += lambda * (cN * grad_N + cS * grad_S + cE * grad_E +
												   cW * grad_W + cF * grad_F + cB * grad_B);

				}

			}
		}
	}

	return denoised;
}


double*** aniso_diffusion_2d(double*** img, float lambda, float kappa, int niter,
						 int nx, int ny, int nz){

	double ***denoised;
	double grad_N, grad_S, grad_E, grad_W;
	double cN, cS, cE, cW, Ksq;
	int jp, jm, ip, im;
	int i, j, k;

	denoised = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(denoised[i][j],img[i][j],sizeof(double)*nz);
	}

	Ksq = pow(kappa, 2);

	while(niter--){

		for(k=0; k<nx; k++){

			for(j=0; j<ny; j++){

				jp = (j==ny-1)? ny-1: j+1;
				jm = (j==0) ? 0 : j-1;

				for(i=0; i<nz; i++){

					ip = (i==nx-1)? nx-1: i+1;
					im = (i==0) ? 0 : i-1;

					grad_N = denoised[k][jm][i] - denoised[k][j][i];
					grad_S = denoised[k][jp][i] - denoised[k][j][i];
					grad_E = denoised[k][j][ip] - denoised[k][j][i];
					grad_W = denoised[k][j][im] - denoised[k][j][i];

					cN = 1.0 / (1.0 + pow(grad_N, 2) / Ksq);
					cS = 1.0 / (1.0 + pow(grad_S, 2) / Ksq);
					cE = 1.0 / (1.0 + pow(grad_E, 2) / Ksq);
					cW = 1.0 / (1.0 + pow(grad_W, 2) / Ksq);
/*
					cN = exp( -pow(grad_N, 2) / Ksq );
					cS = exp( -pow(grad_S, 2) / Ksq );
					cE = exp( -pow(grad_E, 2) / Ksq );
					cW = exp( -pow(grad_W, 2) / Ksq );
					cF = exp( -pow(grad_F, 2) / Ksq );
					cB = exp( -pow(grad_B, 2) / Ksq );
*/
					denoised[k][j][i] += lambda * (cN * grad_N + cS * grad_S + cE * grad_E + cW * grad_W);

				}

			}
		}
	}

	return denoised;
}


double*** show_grad_mag(double*** img, int nx, int ny, int nz){

	double ***grad_mag;
	int i, j, k;

	grad_mag = malloc_image(nx, ny, nz);

	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			memmove(grad_mag[i][j],img[i][j],sizeof(double)*nz);
	}

	for(i=1; i<nx-1; i++){
		for(j=1; j<ny-1; j++){
			for(k=1; k<nz-1; k++)
				grad_mag[i][j][k] = sqrt(pow(img[i][j][k]-img[i-1][j][k],2)+pow(img[i][j][k]-img[i][j-1][k],2)+
										 pow(img[i][j][k]-img[i][j][k-1],2));
		}
	}

	return grad_mag;
}

double*** malloc_image(int nx, int ny, int nz){

	double*** image;
	int i, j;

	image = (double***) malloc(nx*sizeof(double**));

	for(i=0; i<nx; i++){

		image[i] = (double **) malloc(ny*sizeof(double*));

		for(j=0; j<ny; j++)
			image[i][j] = (double*) malloc(nz*sizeof(double));

	}

	return image;
}

void free_image(double ***img, int nx, int ny){

	int i, j;
	for(i=0; i<nx; i++){
		for(j=0; j<ny; j++)
			free(img[i][j]);
		free(img[i]);
	}

	free(img);
	return;
}
