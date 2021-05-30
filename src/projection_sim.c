#include <stdarg.h>
#include <getopt.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "grid.h"

#define NFFT_PRECISION_DOUBLE
#include "nfft3mp.h"

static void show_help(const char *s){

	printf("\nSyntax: %s [options]\n\n",s);
	printf(
		   "Purpose: Simulate serials of tomograms from 3D density arranged in dimensions nx*ny*nz.\n"
		   "         To be noted, rotation axes is along direction Z, perpendicular to paper. X is horizontal while Y is vertical.\n"
		   "         And rotation angle is between X direction and incident rays.\n"
		   "         Thereby, the size of each projection is np*nz.\n"
		   "\n"
		   "-h, --help                    Display this help message.\n"
		   "\n"
		   "-i, --infile=<file>           File from which to get the 3D density.\n"
		   "    --nx=<nx>                 X dimension.\n"
		   "    --ny=<ny>                 Y dimension.\n"
		   "    --nz=<nz>                 Z dimension.\n"
		   "-p, --pad=<np>                Dimension of projection in XY plane.\n"
		   "    --max=<maximum gama>.     Maximum angle of tomograms for simulation.\n"
		   "    --min=<minimum gama>.     Minimum angle of tomograms for simulation.\n"
		   "-t, --interval=<delta angle>. Angular sampling interval.\n"
		   "    --noise.                  Add noise to tomograms.\n"
		   "    --snr=<value>.            Signal to noise level.\n"
		   "-m  --method=<2D1D/3D/DFT>.   Use 2D NFFT + 1D FFT (default) or 3D NFFT (slow) or Discrete Fourier transform.\n"
		   "\n"
		   );
}

typedef enum
{
	FAST,
	SLOW,
	DFT
} ProjectionMethod;

static void add_noise(fftw_complex *slice, int np, int nz, double snr, gsl_rng *rng){

	int i;
	double sum2, sum, sigma;

	sum2=0;
	sum=0;
	//first calculate variance of slice
	for(i=0; i<np*nz; i++){
		sum += creal(slice[i]);
		sum2 += pow(creal(slice[i]), 2);
	}

	sigma = sqrt( sum2/(np*nz) - pow(sum/(np*nz), 2) );

	//add gaussian noise
	for(i=0; i<np*nz; i++){
		slice[i] += gsl_ran_gaussian(rng, sigma/snr);
	}

	return;

}



int main(int argc, char *argv[]){

	int c;
	char *infile = NULL;
	char *method = NULL;
	ProjectionMethod proj_m;
	char *rval;
	double interval = -1;
	double gam_max = 100;
	double gam_min = 100;
	int nx = -1;
	int ny = -1;
	int nz = -1;
	int np = -1;
	int config_noise = 0;
	double snr = -1;

	int nslices;

	FILE *fp;
	fftw_complex *data, *slices;

	gsl_rng *rng;

	/* Long options */
	const struct option longopts[] = {
		{"help",       0,   NULL,        'h'},
		{"infile",     1,   NULL,        'i'},
		{"pad",        1,   NULL,        'p'},
		{"interval",   1,   NULL,        't'},
		{"method",     1,   NULL,        'm'},
		{"nx",         1,   NULL,          2},
		{"ny",         1,   NULL,          3},
		{"nz",         1,   NULL,          4},
		{"max",        1,   NULL,          5},
		{"min",        1,   NULL,          6},
		{"snr",        1,   NULL,          7},
		{"noise",      0,   NULL,          8},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ( ( c = getopt_long(argc, argv, "hi:p:t:m:",longopts,NULL)) != -1 ){

		switch(c) {

		     case 'h':
				 show_help(argv[0]);
				 return 0;

			 case 'i':
				 infile = strdup(optarg);
				 break;

			 case 'p':
				 np = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid X dimension!\n");
					 return 1;
				 }
				 break;

			 case 't':
				 interval = strtod(optarg, &rval);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid angular sampling interval!\n");
					 return 1;
				 }
				 break;

			 case 'm':
				 method = strdup(optarg);
				 break;

			 case 2:
				 nx = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid X dimension!\n");
					 return 1;
				 }
				 break;

			 case 3:
				 ny = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid Y dimension!\n");
					 return 1;
				 }
				 break;

			 case 4:
				 nz = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid Z dimension!\n");
					 return 1;
				 }
				 break;

			 case 5:
				 gam_max = strtod(optarg, &rval);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid maximum collecting angle!\n");
					 return 1;
				 }
				 break;

			 case 6:
				 gam_min = strtod(optarg, &rval);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid minimum collecting angle!\n");
					 return 1;
				 }
				 break;

			 case 7:
				 snr = strtod(optarg, &rval);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid signal to noise ratio!\n");
					 return 1;
				 }
				 break;

			 case 8:
				 config_noise = 1;
				 break;

			 default:
//				 fprintf(stderr,"Unhandled option '%c' \n", c);
				 return 1;
		}
	}

	if( infile == NULL ){
		fprintf(stderr,"You need to provide 3D density with --infile.\n");
		return 1;
	}

	if( nx<0 || ny<0 || nz<0){
		fprintf(stderr, "You need to provide valid dimensions with --nx, --ny, --nz.\n");
		return 1;
	}

	if( np < 0){
		fprintf(stderr, "You need to provide valid projection dimension in XY plane with --pad.\n");
		return 1;
	}

	if( interval < 0){
		fprintf(stderr, "You need to provide valid angular sampling interval with --interval.\n");
		return 1;
	}

	if( gam_max<-90 || gam_max>90 || gam_min<-90 || gam_max>90){
		fprintf(stderr, "You need to provide valid scope of tomogram angles (-90,90] with --max, --min.\n");
		return 1;
	}

	if( method == NULL ){
		fprintf(stderr, "You didn't specify simulating method, so I'm  using"
				" the 2D1D method, which is fastest.\n");
		proj_m = FAST;
	}else if( strcmp(method, "2D1D") == 0 ){
		proj_m = FAST;
		free(method);
	}else if( strcmp(method, "3D") == 0){
		proj_m = SLOW;
		free(method);
	}else if( strcmp(method, "DFT") == 0){
		proj_m = DFT;
		free(method);
	}

	if( config_noise && snr<0){
		fprintf(stderr, "You must provide valid signal to noise ratio with --snr or delete --noise.\n");
		return 1;
	}

	/*read in 3D density*/

	nslices = (int) ( (gam_max - gam_min)/interval + 0.5 ) + 1;

	data = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*nx*ny*nz);
	slices = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*np*nz*nslices);

	/* C-routine array: since higher dimension has faster change, 1D data is actually arranged in data[nz][ny][nx] */
	if(fp = fopen(infile,"r")){

		int i, j, k;
		double real;

		for(k=0; k<nz; k++)
		for(j=0; j<ny; j++)
		for(i=0; i<nx; i++){
			fscanf(fp, "%le\n", &real);
			data[k*nx*ny+j*nx+i] = real;
		}

	}else{
		fprintf(stderr,"Error : Data file %s doesn't exist!\n ", infile);
		return 1;
	}

	fclose(fp);
	free(infile);

	if(config_noise){
		rng = gsl_rng_alloc(gsl_rng_default);
	}


	if(proj_m == FAST){

	/*perform 1D FFT + 2D NFFT*/

		nfft_plan p;
		fftw_plan plan;
		int d = 2;
		int N[2], n[2];
		int M, m;
		double *x;
		int i, j, k, ntotal;

		//construct equal-angle and non-uniform sampling grid. Note that 2D x is arranged in Y --> X
		x = equal_angle_sample_2D(gam_min, gam_max, interval, np, &ntotal);

		//prepare 2D nfft
		M = ntotal;
		N[0] = ny;
		N[1] = nx;
		n[0] = (int) nfft_next_power_of_2(N[0])*2;
		n[1] = (int) nfft_next_power_of_2(N[1])*2;

		nfft_init_guru(&p,2,N,M,n,6,
					   PRE_PHI_HUT | PRE_FULL_PSI | MALLOC_X | MALLOC_F_HAT |
					   MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					   FFTW_MEASURE | FFTW_DESTROY_INPUT);

		//initialize non-uniformly sampled p.x
		for(i=0; i<p.M_total*p.d; i++){
			p.x[i] = x[i];
		}

		free(x);

		if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);

		//1D FFT
		plan = fftw_plan_many_dft(1, &nz, nx*ny, data, NULL, nx*ny, 1,
								  data, NULL, nx*ny, 1, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		//2D NFFT for each slice along z direction
		for(j=0; j<nz; j++){

			//intialize uniformly sampled p.h_hat
			for(i=0; i<nx*ny; i++){

				/* fftshift after 1D fft to move zero frequency to center (because the following nfft2 in frequency is zero-centric)
				 * Also note that the order of storing 1D p.f_hat obeys C routine: higher dimension is first stored.
				 * That is to say, X is first stored, followed by Y, namely p.f_hat[ny][nx]. Also 1D data is arranged in data[nz][ny][nx] */

				int col = i%nx; //along X
				int row = i/nx; //along Y

				p.f_hat[row*nx+col] = data[(j*nx*ny+i+nx*ny*nz/2)%(nx*ny*nz)];

			}

			//conduct 2D nfft to generate equal-angle and uniform Fourier slices
			nfft_trafo(&p);

			//save fourier slices, which is arrayed in slices[nslices][nz][np]
			for(i=0; i<p.M_total; i++){
				int nslc = i/np;
				int ncol = i%np;
				slices[nslc*np*nz+j*np+ncol] = p.f[i];
			}

		}

		nfft_finalize(&p);

		fftw_complex *slice;
		slice = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*np*nz);

		//2D ifft for each slice
		for(j=0; j<nslices; j++){

			//extract one slice
			for(i=0; i<np*nz; i++){
				//ifftshift after nfft2 for moving zero-frequency to corner
				int row = (i/np+nz/2) % nz;
				int col = (i%np+np/2) % np;
				slice[row*np+col] = slices[j*np*nz+i];
			}

			/* fftw for multi-dimension array: if array is arranged in data[nx][ny], then the order for fftw is {nx, ny}.
			 * Here, slice is arranged in slice[nz][np] */

			plan = fftw_plan_dft_2d(nz, np, slice, slice, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			/* fftshift unidimension slice after ifft2 */
			for(i=0; i<np*nz; i++){

				int row = i/np;
				int col = i%np;
				int col_shift = (col+np/2) % np;

				if(col >= np/2) continue;

				//exchange corresponding locations before and after shift
				complex tmp = slice[row*np+col];
				slice[row*np+col] = slice[row*np+col_shift];
				slice[row*np+col_shift] = tmp;

			}

			if(config_noise){
				add_noise(slice, np, nz, snr, rng);
			}

			char filename[1024];

			snprintf(filename,1023,"./files/prj-%i.txt",j);

			fp = fopen(filename, "w");
			for(i=0; i<np*nz; i++)
				fprintf(fp,"%f\n", creal(slice[i])/(np*nz));
			fclose(fp);

		}

		nfft_free(slice);

	}else if(proj_m == SLOW){

	/*perform 3D NFFT*/

		nfft_plan p;
		fftw_plan plan;
		int d = 3;
		int N[3], n[3];
		int M, m;
		double *x;
		int i, j, k, ntotal;

		/* construct equal-angle and non-uniform sampling grid */
		x = equal_angle_sample_3D(gam_min, gam_max, interval, np, nz, &ntotal);

		//prepare 3D nfft
		M = ntotal;
		N[0] = nx;
		N[1] = ny;
		N[2] = nz;
		n[0] = (int) nfft_next_power_of_2(N[0])*2;
		n[1] = (int) nfft_next_power_of_2(N[1])*2;
		n[2] = (int) nfft_next_power_of_2(N[2])*2;

		nfft_init_guru(&p,3,N,M,n,6,
					   PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT |
					   MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					   FFTW_MEASURE | FFTW_DESTROY_INPUT);

		//initialize non-uniformly sampled p.x
		for(i=0; i<p.M_total*p.d; i++){
			p.x[i] = x[i];
		}

		free(x);

		if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);

		/* intialize uniformly sampled p.h_hat
		 *
		 * p.f_hat must be arrayed in [nx]-[ny]-[nz], which means Z is first stored, followed by Y and X. */
		for(i=0; i<nx*ny*nz; i++){

			int col = (i%(nx*ny))%nx; //along X
			int row = (i%(nx*ny))/nx; //along Y
			int slc = i/(nx*ny); //along Z

			p.f_hat[col*ny*nz+row*nz+slc] = data[i];

		}

		/* conduct 3D nfft to generate equal-angle and uniform Fourier slices */
		nfft_trafo(&p);

		/* 2D IFFT for each slice */
		fftw_complex *slice;
		slice = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*np*nz);

		for(j=0; j<nslices; j++){

			for(i=0; i<np*nz; i++){
				//ifftshift after nfft for moving zero-frequency to corner
				int row = (i/np+nz/2) % nz;
				int col = (i%np+np/2) % np;
				slice[row*np+col] = p.f[j*np*nz+i];
			}

			/* fftw for multi-dimension array: if array is arranged in data[nx][ny], then the order for fftw is {nx, ny}.
			 * Here, slice is arranged in slice[nz][np] */

			plan = fftw_plan_dft_2d(nz, np, slice, slice, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			/* fftshift slice after ifft2 */
			for(i=0; i<np*nz/2; i++){

				int row = i/np;
				int col = i%np;
				int row_shift = (row+nz/2) % nz;
				int col_shift = (col+np/2) % np;

				//exchange corresponding locations before and after shift
				complex tmp = slice[row*np+col];
				slice[row*np+col] = slice[row_shift*np+col_shift];
				slice[row_shift*np+col_shift] = tmp;

			}

			if(config_noise){
				add_noise(slice, np, nz, snr, rng);
			}

			char filename[1024];

			snprintf(filename,1023,"./files/prj-%i.txt",j);

			fp = fopen(filename, "w");
			for(i=0; i<np*nz; i++)
				fprintf(fp,"%f\n", creal(slice[i])/(np*nz));
			fclose(fp);

		}

		nfft_finalize(&p);
		nfft_free(slice);

	}else if(proj_m == DFT){

		fftw_plan plan;
		double *x;
		int i, j, k, ntotal;

		/* construct equal-angle and non-uniform sampling grid */
		x = equal_angle_sample_3D(gam_min, gam_max, interval, np, nz, &ntotal);

		fftw_complex *slice;
		slice = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*np*nz);

		for(k=0; k<nslices; k++){

			/*calculate discrete transformation for each slice*/
			for(j=0; j<np*nz; j++){

				slice[j] = 0.0;

				int idx = 3*(k*np*nz+j);

				for(i=0; i<nx*ny*nz; i++){

					int col = (i%(nx*ny))%nx - nx/2; //along X, [-nx/2, nx/2-1]
					int row = (i%(nx*ny))/nx - ny/2; //along Y, [-ny/2, ny/2-1]
					int slc = i/(nx*ny) - nz/2;      //along Z, [-nz/2, nz/2-1]

					complex exponent = cos(2*PI*(col*x[idx+0]+row*x[idx+1]+slc*x[idx+2])) -
						               sin(2*PI*(col*x[idx+0]+row*x[idx+1]+slc*x[idx+2]))*_Complex_I;

					slice[j] += data[i] * exponent;

				}

			}

			/* 2D IFFT for each slice */

			/*ifftshift after nfft for moving zero-frequency to corner*/
			for(i=0; i<np*nz/2; i++){

				int row = i/np;
				int col = i%np;
				int row_shift = (row+nz/2) % nz;
				int col_shift = (col+np/2) % np;

				//exchange corresponding locations before and after shift
				complex tmp = slice[row*np+col];
				slice[row*np+col] = slice[row_shift*np+col_shift];
				slice[row_shift*np+col_shift] = tmp;

			}
			/* fftw for multi-dimension array: if array is arranged in data[nx][ny], then the order for fftw is {nx, ny}.
			 * Here, slice is arranged in slice[nz][np] */

			plan = fftw_plan_dft_2d(nz, np, slice, slice, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			/* fftshift slice after ifft2 */
			for(i=0; i<np*nz/2; i++){

				int row = i/np;
				int col = i%np;
				int row_shift = (row+nz/2) % nz;
				int col_shift = (col+np/2) % np;

				//exchange corresponding locations before and after shift
				complex tmp = slice[row*np+col];
				slice[row*np+col] = slice[row_shift*np+col_shift];
				slice[row_shift*np+col_shift] = tmp;

			}

			if(config_noise){
				add_noise(slice, np, nz, snr, rng);
			}

			char filename[1024];

			snprintf(filename,1023,"./files/prj-%i.txt",k);

			fp = fopen(filename, "w");
			for(i=0; i<np*nz; i++)
				fprintf(fp,"%f\n", creal(slice[i])/(np*nz));
			fclose(fp);

		}


		free(x);
		nfft_free(slice);

	}

	if(config_noise) gsl_rng_free(rng);

	nfft_free(data);
	nfft_free(slices);

	return 0;
}
