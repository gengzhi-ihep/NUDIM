#include "mrcfile.h"

void main(){

	int i, j, k;

	MrcHeader *inhead;

	float *slcdata, *datatmp;

	FILE *fin, *fout;

	char filename[1024];

	inhead = (MrcHeader *) malloc(sizeof(MrcHeader));

	fin = fopen("BBa.ali","r");

	mrc_read_head(fin, inhead);

	printf("Data Dimension: %d %d %d\n", inhead->nx, inhead->ny, inhead->nz);
	printf("Data type (Mode): %d\n", inhead->mode);
	printf("Nstart: %d %d %d\n", inhead->nxstart, inhead->nystart, inhead->nzstart);
	printf("Grid size: %d %d %d\n", inhead->mx, inhead->my, inhead->mz);
	printf("Length: %f %f %f\n", inhead->xlen, inhead->ylen, inhead->zlen);
	printf("Angles: %f %f %f\n", inhead->alpha, inhead->beta, inhead->gamma);
	printf("Map column: %d %d %d\n", inhead->mapc, inhead->mapr, inhead->maps);
	printf("Density: %f %f %f\n", inhead->amin, inhead->amax, inhead->amean);
	printf("Image type: %d", inhead->ispg);
	printf("Space Group: %d\n", inhead->nsymbt);

	slcdata = (float *) malloc(sizeof(float)*inhead->nx*inhead->ny*inhead->nz);
	datatmp = (float *) malloc(sizeof(float)*inhead->nx*inhead->ny);

	for(k=0; k<inhead->nz; k++){
		mrc_read_slice(fin, inhead, k, 'z', &slcdata[k*inhead->nx*inhead->ny]);
	}

	fclose(fin);

	for(k=1; k<inhead->nz; k++){

                double max=-100000, min=100000;

		//shift image by 61 along X (fast axis); 1 along Y (only within 61 columns)
		for(j=0; j<inhead->ny; j++){
			for(i=0; i<inhead->nx; i++){

                                if(slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i] > max)
                                     max = slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i];

                                if(slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i] < min)
                                     min = slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i];

				int col = (i-61+inhead->nx)%inhead->nx;
				int row = (j-1+inhead->ny)%inhead->ny;
				if(i<61){
					datatmp[row*inhead->nx+col] = slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i];
				}else{
					datatmp[j*inhead->nx+col] = slcdata[k*inhead->nx*inhead->ny+j*inhead->nx+i];
				}
                               
			}
		}

		//output
		char filename[1024];
		snprintf(filename,1023,"./files/prj-%i.txt",inhead->nz-1-k);

                printf("%f %f \n", max, min);

		fout = fopen(filename,"w");

		for(i=0; i<inhead->nx*inhead->ny; i++){
			fprintf(fout, "%f\n", max-datatmp[i]+min);
		}

		fclose(fout);

	}

	free(slcdata);
	free(datatmp);


	free(inhead);

	return;

}
