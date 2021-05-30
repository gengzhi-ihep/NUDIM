* NUDIM: a Non-Uniform fast Fourier transform based Dual-space constraint Iterative reconstruction Method

* Author: Dr.Geng Zhi 

* Email address: gengz@ihep.ac.cn

* NUDIM is used for restoring missing-wedge information from limited number of projections and is capable of improving reconstruction resolution for biological electron tomography.

* Please cite the following publication if you use the codes in your work.

* Geng,Z.,She,Z.,Zhou,Q.,Dong,Z.,Zhan,F.,Zhang,H.,Xu,J.H.,Gao,Z.Q. and Dong,Y.H. NUDIM:A non-uniform fast Fourier transform based dual-space constraint iterative reconstruction method 

* in biological electron tomography.

################################
* Prerequisite:

* gcc

* openmpi

* gsl

* nfft

* fftw

################################


################################
* Installation:

* 1. Install gsl, openmpi, nfft and fftw3. Note that all these libraries can be easily installed on Centos 7/ 64-bit oprating system with the yum command as follows:

*    yum install openmpi.x86_64 openmpi-devel.x86_64 gsl.x86_64 gsl-devel.x86_64 fftw.x86_64 fftw-devel.x86_64 nfft.x86_64 nfft-devel.x86_64

* 2. Modify the file 'Makefile' to assign your library path and header-file path to openmpi, gsl, nfft and fftw.

* 3. Then run 'make all'. This will creat three executable files in directory '/bin', including 'projection_sim', 'reconstruction', 'reconstruction_parallel'.

################################

################################
* Program description:

1. projection_sim:

   Syntax: ./projection_sim [options]

   Purpose: Simulate serials of tomograms from 3D density arranged in dimensions nx*ny*nz.
            To be noted, rotation axes is along direction Z, perpendicular to paper. X is horizontal while Y is vertical.
            And rotation angle is between X direction and incident rays.
            Thereby, the size of each projection is np*nz.

   -h, --help                    Display this help message.

   -i, --infile=<file>           File from which to get the 3D density.
       --nx=<nx>                 X dimension.
       --ny=<ny>                 Y dimension.
       --nz=<nz>                 Z dimension.
   -p, --pad=<np>                Dimension of projection in XY plane.
       --max=<maximum gama>.     Maximum angle of tomograms for simulation.
       --min=<minimum gama>.     Minimum angle of tomograms for simulation.
   -t, --interval=<delta angle>. Angular sampling interval.
       --noise.                  Add noise to tomograms.
       --snr=<value>.            Signal to noise level.
   -m  --method=<2D1D/3D/DFT>.   Use 2D NFFT + 1D FFT (default) or 3D NFFT (slow) or Discrete Fourier transform.

2. reconstruction:

   Syntax: ./reconstruction [options]

   Purpose: Reconstruct 3D object from series of equally-angled tomograms.
            To be noted, rotation axes is along direction Z, perpendicular to paper. X is horizontal while Y is vertical.
            And rotation angle is between X direction and incident rays.
            Thereby, the size of each projection is nx*nz , and nx is first ordered, namely projection[nz][nx].

   -h, --help                    Display this help message.

       --prefix=<prefix name>.   Prefix name of series of projections.
       --nx=<nx>                 X dimension size of projection.
       --nz=<nz>                 Z dimension size of projection.
   -p, --pad=<np>                Padded size along radial direction, which must be at least larger than nx to achieve oversampling.
       --max=<maximum gama>.     Maximum angle of tomograms.
       --min=<minimum gama>.     Minimum angle of tomograms.
   -t, --interval=<delta angle>. Angular sampling interval.
       --padangle=<pad interval>.Filling unknown region every padangle, which should be divided by the above interval.
       --weight.                 Add weight to fourier slices to compensate for uneven sampling density among different frequencies.
       --denoise.                Filtering 3D density at each iteration.
       --method=<DIFFUSION/TV>.  Use anisotropic diffusion or total variation to denoise 3D density.
       --nhio=<nhio>.            Number of HIO iterations.
       --ner=<ner>.              Number of ER iterations.
       --height.                 Sample height for constraint.
       --random.                 Padding unknown projections with randoms or zeros.
       --print.                  Print out R factor with each cycle of iteration.


################################
* How to use:

* Several example scripts are provided in /scripts.

* There are also some example datasets in /examples for test.

################################
