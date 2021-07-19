NUDIM:A Non-Uniform fast Fourier transform based Dual-space constraint Iterative reconstruction Method
===================
Author: Zhi Geng | email: gengz@ihep.ac.cn

Description:
-------------------

`NUDIM` is a software developed to restore missing-wedge information from limited number of projections in biological electron tomography, which has the following features:      
* It is based on dual-space constraint iteration which iteratively imposes constraints between real space and Fourier space.  
* It makes use of *NFFT* to transform between real space and reciprocal space, thus circumventing interpolation induced errors.  
* Totol variation based image denoising is further incorporated into dual-space constraint to suppress noise.  

Please cite the following publication if you use ***NUDIM*** in your work: 😊
> Geng,Z.,She,Z.,Zhou,Q.,Dong,Z.,Zhan,F.,Zhang,H.,Xu,J.H.,Gao,Z.Q. & Dong,Y.H. NUDIM:A non-uniform fast Fourier transform based dual-space constraint iterative reconstruction method in biological electron tomography. Journal of Structural Biology

Requirements:
------------------
Before compiling source codes, the following binary libraries are required:

* gcc
* openmpi
* gsl
* nfft
* fftw

Installation:
------------------

1. Install library: `gsl`, `openmpi`, `nfft` and `fftw3`.  

 * Note that all these libraries can be easily installed through `yum` on our test operating system: *Centos7-64 bits*:  
 
```
yum install openmpi.x86_64 openmpi-devel.x86_64 gsl.x86_64 gsl-devel.x86_64 fftw.x86_64 fftw-devel.x86_64 nfft.x86_64 nfft-devel.x86_64
```

2. Modify `Makefile` to specify library path and include path to `openmpi`, `gsl`, `nfft` and `fftw`.  

3. Run `make all`.
  
 * The following 3 binary files will be created in directory `bin`:
   
   `projection_sim`,   `reconstruction`,   `reconstruction_parallel`  

4. Finish. Now you can enjoy your test! 😃

Usage:
-------------------
1. Simulation of tilt series from 3D density map in text format.  

* Example: using the supplied 4v6i.txt map, projection_sim can be run with:  

```
projection_sim -i 4v6i.txt -p 128 --nx=128 --ny=128 --nz=128 --max=60 --min=-60 -t 3 -m 2D1D --noise --snr=1
```

`Purpose`: Simulate tilt series from 3D density map with dimensions of (nx,ny,nz). Rotation axis is along Z axis perpendicular to paper, X axis is horizontal. Each projection is with a size of (np,nz). np can be oversampled to make a boundary.  

`Note`: You must create a directory `files` under current path to save simluated projections before running.

* Options you may set are:

```
   -h, --help                    Display this help message.

   -i, --infile=<file>           File from which to get the 3D density.
       --nx=<nx>                 X dimension size.
       --ny=<ny>                 Y dimension size.
       --nz=<nz>                 Z dimension size.
   -p, --pad=<np>                Size of projection in XY slice perpendicular to tilt axis.
       --max=<maximum gama>.     Maximum angle of tilt series for simulation.
       --min=<minimum gama>.     Minimum angle of tilt series for simulation.
   -t, --interval=<delta angle>. Angular sampling interval.
       --noise.                  Add noise to tilt series.
       --snr=<value>.            Signal to noise level.
   -m  --method=<2D1D/3D/DFT>.   Use 2D NFFT + 1D FFT (default and very fast) or 3D NFFT (medium speed) or Discrete Fourier transform(very slow, most accurate).
```

2. 3D reconstruction with ***NUDIM***.

* Examples: provided projections with prefix of 'prj' saved in `files`, ***NUDIM*** can be run with:

```
reconstruction --prefix=prj --pad=128 --nx=128 --nz=128 --max=60 --min=-60 -t 3 --padangle=1 --nhio=0 --ner=20 --denoise --method=TV --print --random --weight
```

`Purpose`: Reconstruct 3D object from tilt series of equally-angled projections.

* Options you may set are:

```
   -h, --help                    Display this help message.

       --prefix=<prefix name>.   Prefix name of projections saved in files.
       --nx=<nx>                 X dimension size of projection.
       --nz=<nz>                 Z dimension size of projection.
   -p, --pad=<np>                Padded size along radial direction, which can be larger than nx to achieve oversampling.
       --max=<maximum gama>.     Maximum angle of tilt series.
       --min=<minimum gama>.     Minimum angle of tilt series.
   -t, --interval=<delta angle>. Experimental angular sampling interval.
       --padangle=<pad interval>.Filling unknown region with an angular interval of padangle, which must be divisible by the above experimental angular interval.
       --weight.                 Add weight to fourier slices to compensate for uneven sampling density in different frequencies.
       --denoise.                Filtering/denoising 3D density at each iteration.
       --method=<DIFFUSION/TV>.  Use anisotropic diffusion (DIFFUSION) or total variation minimization (TV) to denoise 3D density.
       --nhio=<nhio>.            Number of HIO iterations.
       --ner=<ner>.              Number of ER iterations.
       --height=<nheight>.       Sample height for constraint (in pixel).
       --center=<ncenter>.       Sample center for constraint (in pixel).
       --random.                 Padding unknown projections with random values or zeros during initilization of constructed Fourier grid.
       --print.                  Print out R factor with each cycle of iteration.
```

3. 3D reconstruction with ***NUDIM*** with Parallel Computing.

* Examples: run the same example as above with 4 CPU kernels:

```
/usr/lib64/openmpi/bin/mpirun -n 5 reconstruction_parallel --prefix=prj --pad=128 --nx=128 --nz=128 --max=60 --min=-60 -t 3 --padangle=1 --nhio=0 --ner=20 --denoise --method=TV --print --random --weight
```
