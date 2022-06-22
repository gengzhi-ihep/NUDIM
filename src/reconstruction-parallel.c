#include <stdarg.h>
#include <getopt.h>
#include <complex.h>
#include <gsl/gsl_rng.h>

#include "grid.h"
#include "denoise.h"

#define NFFT_PRECISION_DOUBLE
#include "nfft3mp.h"

#include "mpi.h"

static void show_help(const char *s){

	printf("\nSyntax: %s [options]\n\n",s);
	printf(
		   "Purpose: Reconstruct 3D object from series of equally-angled tomograms.\n"
		   "         To be noted, rotation axes is along direction Z, perpendicular to paper. X is horizontal while Y is vertical.\n"
		   "         And rotation angle is between X direction and incident rays.\n"
		   "         Thereby, the size of each projection is nx*nz , and nx is first ordered, namely projection[nz][nx].\n"
		   "\n"
		   "-h, --help                    Display this help message.\n"
		   "\n"
		   "    --prefix=<prefix name>.   Prefix name of series of projections.\n"
		   "    --nx=<nx>                 X dimension size of projection.\n"
		   "    --nz=<nz>                 Z dimension size of projection.\n"
		   "-p, --pad=<np>                Padded size along radial direction, which must be at least larger than nx to achieve oversampling.\n"
		   "    --max=<maximum gama>.     Maximum angle of tomograms.\n"
		   "    --min=<minimum gama>.     Minimum angle of tomograms.\n"
		   "-t, --interval=<delta angle>. Angular sampling interval.\n"
		   "    --padangle=<pad interval>.Filling unknown region every padangle, which should be divided by the above interval.\n"
                   "    --weight.                 Add weight to fourier slices to compensate for uneven sampling density among different frequencies.\n"
		   "    --denoise.                Filtering 3D density at each iteration.\n"
		   "    --method=<DIFFUSION/TV>.  Use anisotropic diffusion or total variation to denoise 3D density.\n"
		   "    --nhio=<nhio>.            Number of HIO iterations.\n"
		   "    --ner=<ner>.              Number of ER iterations.\n"
		   "    --height.                 Sample height for constraint.\n"
		   "    --random.                 Padding unknown projections with randoms or zeros.\n"
		   "    --print.                  Print out R factor with each cycle of iteration.\n"
		   "\n"
		  );
}

typedef enum
{
	DIFFUSION,
	TV
} DenoiseMethod;

int main(int argc, char *argv[]){

	int c;
	char *prefix = NULL;
	char *method = NULL;
	DenoiseMethod denoise_m;
	char *rval;
	double interval = -1;
	double padangle = -1;
	double gam_max = 100;
	double gam_min = 100;
	int nx = -1;
	int nz = -1;
        int nheight = -1;
	int np = -1;
	int nhio = -1;
	int ner = -1;
	int config_denoise = 0;
	int config_random = 0;
	int config_weight = 0;
	int config_print = 0;
	int config_height = 0;

	int nslices;

	fftw_complex *subslices, *data;
	double *tomograms_real, *tomograms_imag;
	double *subslices_real, *subslices_imag;
	double *subdata_real, *subdata_imag;
	double *alldata_real, *alldata_imag;
	double *subimg;

	double ***img, ***img_last, ***denoised, ***grad_mag;
	int *mask;

	nfft_plan p;
	nfft_plan p1;
	solver_plan_complex ip;
	fftw_plan plan;
	int d = 2;
	int N[2], n[2];
	int M, m;
	double *x;
	int i, j, k, nn, ntotal, ntotal1;
	double diff1, diff2, tot_diff1, tot_diff2;
	double start, end;

	int rank, size;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	start = MPI_Wtime();

	/* Long options */
	const struct option longopts[] = {
		{"help",       0,   NULL,          'h'},
		{"pad",        1,   NULL,          'p'},
		{"interval",   1,   NULL,          't'},
		{"prefix",     1,   NULL,            2},
		{"nx",         1,   NULL,            3},
		{"nz",         1,   NULL,            4},
		{"max",        1,   NULL,            5},
		{"min",        1,   NULL,            6},
		{"method",     1,   NULL,            7},
		{"nhio",       1,   NULL,            8},
		{"ner",        1,   NULL,            9},
		{"denoise",    0,   NULL,           10},
		{"random",     0,   NULL,           11},
		{"print",      0,   NULL,           12},
		{"padangle",   1,   NULL,           13},
		{"height",     1,   NULL,           14},
		{"weight",     0,   NULL,           15},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ( ( c = getopt_long(argc, argv, "hp:t:",longopts,NULL)) != -1 ){

		switch(c) {

		     case 'h':
				 show_help(argv[0]);
				 return 0;

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

			 case 2:
				 prefix = strdup(optarg);
				 break;

			 case 3:
				 nx = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid X dimension!\n");
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
				 method = strdup(optarg);
				 break;

			 case 8:
				 nhio = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid number of HIO!\n");
					 return 1;
				 }
				 break;

			 case 9:
				 ner = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid number of ER!\n");
					 return 1;
				 }
				 break;

			 case 10:
				 config_denoise = 1;
				 break;

			 case 11:
				 config_random = 1;
				 break;

			 case 12:
				 config_print = 1;
				 break;

			 case 13:
				 padangle = strtod(optarg, &rval);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid padded angular interval!\n");
					 return 1;
				 }
				 break;

			 case 14:
                                 config_height = 1;
				 nheight = strtol(optarg, &rval, 10);
				 if( *rval != '\0'){
					 fprintf(stderr,"Invalid sample height!\n");
					 return 1;
				 }
				 break;

			 case 15:
				 config_weight = 1;
				 break;

			 default:
//				 fprintf(stderr,"Unhandled option '%c' \n", c);
				 return 1;
		}

	}


	if( prefix == NULL ){
		fprintf(stderr,"You need to provide prefix name of tomogram files with --prefix.\n");
		return 1;
	}

	if( nx<0 || nz<0){
		fprintf(stderr, "You need to provide valid projection size with --nx, --nz.\n");
		return 1;
	}

	if( np < 0){
		fprintf(stderr, "You need to provide valid projection dimension in XY plane with --pad.\n");
		return 1;
	}

	if( nhio < 0){
		fprintf(stderr, "You need to provide valid nhio with --nhio.\n");
		return 1;
	}

	if( ner < 0){
		fprintf(stderr, "You need to provide valid ner with --ner.\n");
		return 1;
	}

	if( interval < 0){
		fprintf(stderr, "You need to provide valid angular sampling interval with --interval.\n");
		return 1;
	}

	if( padangle < 0 || interval/padangle-(int)(interval/padangle+0.1)>1e-5){
		fprintf(stderr, "You need to provide valid padded angular sampling interval with --padangle.\n");
		return 1;
	}

	if( gam_max<-90 || gam_max>90 || gam_min<-90 || gam_max>90){
		fprintf(stderr, "You need to provide valid scope of tomogram angles (-90,90] with --max, --min.\n");
		return 1;
	}

	if( config_denoise && method == NULL){
		fprintf(stderr, "You must provide valid denoising method or delete --denoise.\n");
		return 1;
	}

	if( config_height && nheight <= 0){
		fprintf(stderr, "You must provide valid sample height --height.\n");
		return 1;
	}

	if(method == NULL){
	}else if( strcmp(method, "DIFFUSION") == 0 ){
		denoise_m = DIFFUSION;
		free(method);
	}else if( strcmp(method, "TV") == 0 ){
		denoise_m = TV;
		free(method);
	}

	nslices = (int) ( (gam_max - gam_min)/interval + 0.5 ) + 1;

	if(size==1 || size>nz || nz%(size-1)) MPI_Abort(MPI_COMM_WORLD, 1);

	/* x is arranged in sequence of X --> Z */
	x = uniform_sample_2D(np, nz, &ntotal);

	/* tomograms store fourier slice of each projection */
	tomograms_real = (double*) malloc(sizeof(double)*ntotal*nslices);
	tomograms_imag = (double*) malloc(sizeof(double)*ntotal*nslices);

	if(!rank){

		/* prepare 2D NFFT*/
		M = ntotal;
		N[0] = nx;
		N[1] = nz;
		n[0] = (int) nfft_next_power_of_2(N[0])*2;
		n[1] = (int) nfft_next_power_of_2(N[1])*2;

		nfft_init_guru(&p,2,N,M,n,6,
					   PRE_PHI_HUT | PRE_FULL_PSI | MALLOC_X | MALLOC_F_HAT |
					   MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
					   FFTW_MEASURE | FFTW_DESTROY_INPUT);

		/* initialize uniformly sampled p.x */
		for(i=0; i<p.M_total*p.d; i++){
			p.x[i] = x[i];
		}


		if(p.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);

		/* read in tomogram and calculate oversampled fourier slices one by one */
		for(j=0; j<nslices; j++){

			char filename[1024];
			FILE *fp;

			double tmp;

			snprintf(filename, 1023, "./files/%s-%i.txt", prefix, j);

			fp = fopen(filename, "r");

			for(i=0; i<nx*nz; i++){

				fscanf(fp, "%le\n", &tmp);

				/* since p.x is arrayed in X-->Z, p.f_hat must be inverted: Z-->X,
				 * * that is, Z is first stored, followed by X */
				int col = i%nx; //along X
				int row = i/nx; //along Z

				p.f_hat[col*nz+row] = tmp;

			}

			fclose(fp);

			//conduct 2D nfft to generate oversampled fourier slices of each tomogram.
			nfft_trafo(&p);

			//save oversampled fourier slices
			for(i=0; i<p.M_total; i++){
				tomograms_real[j*p.M_total+i] = creal(p.f[i]);
				tomograms_imag[j*p.M_total+i] = cimag(p.f[i]);
			}

		}

		free(prefix);
		nfft_finalize(&p);

	}

	free(x);

	MPI_Bcast(tomograms_real, ntotal*nslices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(tomograms_imag, ntotal*nslices, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double pad_max, pad_min;
	int ntmp;

	ntmp = (int) ( (-90-gam_min)/interval );
	pad_min = gam_min + ntmp*interval;

	ntmp = (int) ( (90-gam_min)/interval );
	pad_max = gam_min + ntmp*interval;

	/* construct full-range and equally-angle non-uniform sampling grid */
	x = equal_angle_sample_2D(pad_min, pad_max, padangle, np, &ntotal1);

	nslices = ntotal1/np;

	/* prepare 2D NFFT for slices */
	M = ntotal1;
	N[0] = np;
	N[1] = np;
	n[0] = (int) nfft_next_power_of_2(N[0])*2;
	n[1] = (int) nfft_next_power_of_2(N[1])*2;

	nfft_init_guru(&p1,2,N,M,n,6,
				   PRE_PHI_HUT | PRE_FULL_PSI | MALLOC_X | MALLOC_F_HAT |
				   MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE,
				   FFTW_MEASURE | FFTW_DESTROY_INPUT);

	/* initialize non-uniformly sampled p1.x */
	for(i=0; i<p1.M_total*p1.d; i++){
		p1.x[i] = x[i];
	}

	free(x);

	if(p1.flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p1);

	unsigned infft_flags = CGNR;

        if(config_weight){
            infft_flags = infft_flags | PRECOMPUTE_WEIGHT;
        }

	solver_init_advanced_complex(&ip, (nfft_mv_plan_complex *) &p1, infft_flags);

        if(config_weight){
            //generate weights in frequency domain: higher frequency have higher weights.
            for(j=0; j<nslices; j++){
                for(i=0; i<np; i++){
                    if(i == np/2){
                          ip.w[j*np+i] = PI/(4*nslices*pow(np/2,2));
                    }else{
                          ip.w[j*np+i] = padangle/180*PI*abs(i-np/2)/pow(np/2,2);
                    }
                }
            }
        }

	/* slices hoard series of Fourier slices along Z direction. */
	subslices_real = (double*) nfft_malloc(sizeof(double)*ntotal1*nz/(size-1));
	subslices_imag = (double*) nfft_malloc(sizeof(double)*ntotal1*nz/(size-1));

	/* mask is used to label known slices */
	mask = (int*) malloc(sizeof(int)*ntotal1*nz/(size-1));

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

	//allocate data for other jobs
	if(!rank){

		/* reorder from series of tomograms to X-Y-Z slices in order to facilitate 2D NFFT(in XY plane) + 1D FFT(along Z) */
		for(k=0; k<size-1; k++){

			for(j=0; j<nz/(size-1); j++){

				for(i=0; i<ntotal1; i++){

					int nslc = i/np;
					double angle = pad_min + nslc*padangle;
					int span = (int) (interval/padangle+0.5);
					int seq = (int) ( (angle-gam_min)/padangle+0.5 );

					if( (angle < gam_min) || (angle > gam_max) || (seq%span != 0) ){
						/* Unknown reciprocal region */
						mask[j*ntotal1+i] = 0;

						if(config_random){
							/*random padding*/
							subslices_real[j*ntotal1+i] = gsl_rng_uniform(rng)*10;
							subslices_imag[j*ntotal1+i] = 0;
						}else{
							/*zero padding*/
							subslices_real[j*ntotal1+i] = 0;
							subslices_imag[j*ntotal1+i] = 0;
						}

                                        }else if( abs((int)angle) <1e-5 ){
                                                /* cross validation */
                                                mask[j*ntotal1+i] = 0;

						if(config_random){
							/*random padding*/
							subslices_real[j*ntotal1+i] = gsl_rng_uniform(rng)*10;
							subslices_imag[j*ntotal1+i] = 0;
						}else{
							/*zero padding*/
							subslices_real[j*ntotal1+i] = 0;
							subslices_imag[j*ntotal1+i] = 0;
						}

					}else{
						/* Known reciprocal region */
						mask[j*ntotal1+i] = 1;

						int order = (int) ( (angle-gam_min)/interval+0.5 ); //slice number in tomograms
						int col = i%np;

						subslices_real[j*ntotal1+i] = tomograms_real[order*ntotal+(k*nz/(size-1)+j)*np+col];
						subslices_imag[j*ntotal1+i] = tomograms_imag[order*ntotal+(k*nz/(size-1)+j)*np+col];

					}

				}


			}

			MPI_Send(subslices_real, ntotal1*nz/(size-1), MPI_DOUBLE, k+1, 99, MPI_COMM_WORLD);
			MPI_Send(subslices_imag, ntotal1*nz/(size-1), MPI_DOUBLE, k+1, 100, MPI_COMM_WORLD);
			MPI_Send(mask, ntotal1*nz/(size-1), MPI_INT, k+1, 101, MPI_COMM_WORLD);

		}

	}else{

		MPI_Recv(subslices_real, ntotal1*nz/(size-1), MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(subslices_imag, ntotal1*nz/(size-1), MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &status);
		MPI_Recv(mask, ntotal1*nz/(size-1), MPI_INT, 0, 101, MPI_COMM_WORLD, &status);

	}

	gsl_rng_free(rng);
	free(tomograms_real);
	free(tomograms_imag);

	MPI_Barrier(MPI_COMM_WORLD);

	/* Calculate initial 3D density using Inverse NFFT */

	//prepare 1D IFFT
	subdata_real = (double*) malloc(sizeof(double)*np*np*nz/(size-1));
	subdata_imag = (double*) malloc(sizeof(double)*np*np*nz/(size-1));

	subimg = (double*) malloc(sizeof(double)*np*np*nz/(size-1));

	if(rank){

		subslices = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*ntotal1*nz/(size-1));

		for(k=0; k<ntotal1*nz/(size-1); k++){
			subslices[k] = subslices_real[k]+subslices_imag[k]*I;
		}

		//calculate 2D IFFT for each slice
		for(j=0; j<nz/(size-1); j++){

			//read in each slice of nfft data
			for(i=0; i<p1.M_total; i++)
				ip.y[i] = subslices[j*ntotal1+i];

			//init some guess
			for(i=0; i<p1.N_total; i++)
				ip.f_hat_iter[i] = 0.0;

			//conduct 2D IFFT for each slice
			solver_before_loop_complex(&ip);
			for(i=0; i<10; i++){
//				printf("Slice: %d Iteration: %d Residual: %f\n", j, i, ip.dot_r_iter);
				solver_loop_one_step_complex(&ip);
			}

			//copy INFFT data for 1D IFFT
			for(i=0; i<p1.N_total; i++){
				subdata_real[j*np*np+i] = creal(ip.f_hat_iter[i]);
				subdata_imag[j*np*np+i] = cimag(ip.f_hat_iter[i]);
			}

		}

		nfft_free(subslices);

	}

	if(!rank){

		alldata_real = (double*) malloc(sizeof(double)*np*np*nz/(size-1)*size);
		alldata_imag = (double*) malloc(sizeof(double)*np*np*nz/(size-1)*size);

	}

	MPI_Gather(subdata_real, np*np*nz/(size-1), MPI_DOUBLE, alldata_real, np*np*nz/(size-1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(subdata_imag, np*np*nz/(size-1), MPI_DOUBLE, alldata_imag, np*np*nz/(size-1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(!rank){

		data = (fftw_complex*) nfft_malloc(sizeof(fftw_complex)*np*np*nz);

		//ifftshift after infft2 for moving zero-frequency to corner(used for ifft and consistent with ifft(ifftshift))
		for(i=0; i<np*np*nz; i++){
			data[(i+np*np*nz/2)%(np*np*nz)] = alldata_real[np*np*nz/(size-1)+i]+alldata_imag[np*np*nz/(size-1)+i]*I;
		}

		//1D IFFT
		plan = fftw_plan_many_dft(1, &nz, np*np, data, NULL, np*np, 1,
								  data, NULL, np*np, 1, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		//save initial map
		img = malloc_image(nz, np, np);
		img_last = malloc_image(nz, np, np);

		//write out recovered data
		FILE *fp = fopen("recovered-0-parallel.txt","w");

                //FIXME:For some reason, img is somehow flipped at first sight. To make it look good, fftshift along Z is used on img.
		for(j=0; j<nz; j++){
		    for(i=0; i<p1.N_total; i++){
			img[(j+nz/2)%nz][i/np][i%np] = creal(data[j*np*np+i])/nz;
			fprintf(fp, "%le\n", creal(data[(j*np*np+i+np*np*nz/2)%(np*np*nz)])/nz);
                    }
                }
		fclose(fp);

	}

	MPI_Barrier(MPI_COMM_WORLD);


	/* HIO/ER Iteration Begins here */

	for(nn=0; nn<nhio+ner; nn++){

		if(!rank){

			//save 3D image into 1D data
			for(j=0; j<nz; j++){
			for(i=0; i<p1.N_total; i++){
				data[j*np*np+i] = img[(j+nz/2)%nz][i/np][i%np]; //FIXME:Since img is flipped, ifftshift is used again to coincide with NFFT.
				img_last[j][i/np][i%np] = img[j][i/np][i%np];
			}
			}

			//perform 1D FFT
			plan = fftw_plan_many_dft(1, &nz, np*np, data, NULL, np*np, 1,
									  data, NULL, np*np, 1, FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			//allocate subjob datas
			for(k=0; k<size-1; k++){

				for(j=0; j<nz/(size-1); j++){

					for(i=0; i<p1.N_total; i++){
						subdata_real[j*np*np+i] = creal(data[((k*nz/(size-1)+j)*np*np+np*np*nz/2+i)%(np*np*nz)]); //fftshift
						subdata_imag[j*np*np+i] = cimag(data[((k*nz/(size-1)+j)*np*np+np*np*nz/2+i)%(np*np*nz)]);

                                                //FIXME:Since img is flipped, ifftshift is used again to coincide with NFFT.
						subimg[j*np*np+i] = img[(k*nz/(size-1)+j+nz/2)%nz][i/np][i%np];
					}
				}

				MPI_Send(subdata_real, np*np*nz/(size-1), MPI_DOUBLE, k+1, 102, MPI_COMM_WORLD);
				MPI_Send(subdata_imag, np*np*nz/(size-1), MPI_DOUBLE, k+1, 103, MPI_COMM_WORLD);
				MPI_Send(subimg, np*np*nz/(size-1), MPI_DOUBLE, k+1, 104, MPI_COMM_WORLD);

			}

		}else{

			MPI_Recv(subdata_real, np*np*nz/(size-1), MPI_DOUBLE, 0, 102, MPI_COMM_WORLD, &status);
			MPI_Recv(subdata_imag, np*np*nz/(size-1), MPI_DOUBLE, 0, 103, MPI_COMM_WORLD, &status);
			MPI_Recv(subimg, np*np*nz/(size-1), MPI_DOUBLE, 0, 104, MPI_COMM_WORLD, &status);

		}

		MPI_Barrier(MPI_COMM_WORLD);

		diff1 = 0;
		diff2 = 0;

		if(rank){

			//Slice by Slice
			for(j=0; j<nz/(size-1); j++){

				for(i=0; i<p1.N_total; i++)
					p1.f_hat[i] = subdata_real[j*np*np+i]+subdata_imag[j*np*np+i]*I;

				//perform NFFT for each slice
				nfft_trafo(&p1);

				//enforce reciprocal constraints
				for(i=0; i<p1.M_total; i++){
					if(mask[j*p1.M_total+i] > 0.5){
						ip.y[i] = subslices_real[j*p1.M_total+i]+subslices_imag[j*p1.M_total+i]*I;
						diff1 += pow(creal(p1.f[i])-subslices_real[j*p1.M_total+i],2)+pow(cimag(p1.f[i])-subslices_imag[j*p1.M_total+i],2);
						diff2 += pow(creal(p1.f[i])+subslices_real[j*p1.M_total+i],2)+pow(cimag(p1.f[i])+subslices_imag[j*p1.M_total+i],2);
					}else{
						ip.y[i] = p1.f[i];
					}
				}

				//perform inverse NFFT

				//first init some guess
				for(i=0; i<p1.N_total; i++)
					ip.f_hat_iter[i] = subimg[j*np*np+i];

				//conduct 2D INFFT for each slice
				solver_before_loop_complex(&ip);
				for(i=0; i<10; i++){
					solver_loop_one_step_complex(&ip);
				}

				//copy INFFT data for 1D IFFT
				for(i=0; i<p1.N_total; i++){
					subdata_real[j*np*np+i] = creal(ip.f_hat_iter[i]);
					subdata_imag[j*np*np+i] = cimag(ip.f_hat_iter[i]);
				}

			}

		}

		MPI_Gather(subdata_real, np*np*nz/(size-1), MPI_DOUBLE, alldata_real, np*np*nz/(size-1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(subdata_imag, np*np*nz/(size-1), MPI_DOUBLE, alldata_imag, np*np*nz/(size-1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Reduce(&diff1, &tot_diff1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&diff2, &tot_diff2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		if(!rank){

			if(config_print)
				printf("Iteration: %d Measured Residual: %f\n", nn, sqrt(tot_diff1/tot_diff2));

			//perform 1D IFFT
			for(i=0; i<np*np*nz; i++){
				data[(i+np*np*nz/2)%(np*np*nz)] = alldata_real[np*np*nz/(size-1)+i]+alldata_imag[np*np*nz/(size-1)+i]*I;
			}
			plan = fftw_plan_many_dft(1, &nz, np*np, data, NULL, np*np, 1,
									  data, NULL, np*np, 1, FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			for(j=0; j<nz; j++){
			for(i=0; i<p1.N_total; i++){

				int yi = i/np;
				int xi = i%np;

                                //FIXME:For some reason, img is somehow flipped at first sight. To make it look good, fftshift along Z is used on img.
				img[j][yi][xi] = creal(data[(j*np*np+i+nz*np*np/2)%(nz*np*np)])/nz;

			}
			}

			//denoising
			if(config_denoise){

				if(denoise_m == DIFFUSION){
					denoised = aniso_diffusion_2d(img, 0.4, 0.005, 80, nz, np, np);
				}else{
					denoised = total_variation_discrete(img, 0.2, 40, nz, np, np);
				}

				//show initial gradient map for visualization of denoising effect
				if(nn == 0){

					FILE *fp;
					grad_mag = show_grad_mag(img, nz, np, np);
					fp = fopen("grad_before.txt", "w");
					for(i=0; i<nz; i++)
					for(j=0; j<np; j++)
					for(k=0; k<np; k++)
						fprintf(fp,"%f\n", grad_mag[i][j][k]);
					fclose(fp);
					free_image(grad_mag, nz, np);

					grad_mag = show_grad_mag(denoised, nz, np, np);
					fp = fopen("grad_after.txt", "w");
					for(i=0; i<nz; i++)
					for(j=0; j<np; j++)
					for(k=0; k<np; k++)
						fprintf(fp,"%f\n", grad_mag[i][j][k]);
					fclose(fp);
					free_image(grad_mag, nz, np);

				}

				for(i=0; i<nz; i++)
				for(j=0; j<np; j++)
				for(k=0; k<np; k++){
					img[i][j][k] = denoised[i][j][k];
				}

				free_image(denoised, nz, np);


			}

			//now it's time to enforce density constraints
			for(j=0; j<nz; j++){
			for(i=0; i<p1.N_total; i++){

				int yi = i/np;
				int xi = i%np;

				if(config_height){

					if(yi>=(np-nheight)/2 && yi<=(np+nheight)/2){
						if(img[j][yi][xi] < 0)
							img[j][yi][xi] = (nn<nhio) ? img_last[j][yi][xi]-0.9*img[j][yi][xi] : 0;
					}else{
						img[j][yi][xi] = (nn<nhio) ? img_last[j][yi][xi]-0.9*img[j][yi][xi] : 0;
					}

				}else{

					if( (yi>=(np-nx)/2) && (yi<(np+nx)/2) &&
					   (xi>=(np-nx)/2) && (xi<(np+nx)/2) ){
						if(img[j][yi][xi] < 0)
							img[j][yi][xi] = (nn<nhio) ? img_last[j][yi][xi]-0.9*img[j][yi][xi] : 0;
					}else{
						img[j][yi][xi] = (nn<nhio) ? img_last[j][yi][xi]-0.9*img[j][yi][xi] : 0;
					}

				}

			}
			}

		}


	}

	if(!rank){

		//save recovered density, neglecting padded regions
		FILE *fp = fopen("recovered-parallel.txt","w");
		for(i=0; i<nz; i++)
		for(j=0; j<np; j++)
		for(k=0; k<np; k++){

			if( j>=(np+nx)/2 || j<(np-nx)/2 ||
			    k>=(np+nx)/2 || k<(np-nx)/2 ) continue;

			fprintf(fp, "%f\n", img[i][j][k]);

		}

                fclose(fp);

	}


	if(!rank){

		free(alldata_real);
		free(alldata_imag);
		nfft_free(data);
		free_image(img, nz, np);
		free_image(img_last, nz, np);

	}

	free(subslices_real);
	free(subslices_imag);
	free(subdata_real);
	free(subdata_imag);
	free(subimg);

	free(mask);

	solver_finalize_complex(&ip);
	nfft_finalize(&p1);

	end = MPI_Wtime();

	if(!rank) printf("Total time used : %f seconds.\n", end-start);

	MPI_Finalize();

	return 0;

}

