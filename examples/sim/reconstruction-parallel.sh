#!/bin/bash
/usr/lib64/openmpi/bin/mpirun -n 9 ../../bin/reconstruction_parallel --prefix=prj --pad=128 --nx=128 --nz=128 \
--max=60 --min=-60 -t 3 --padangle=1 --nhio=0 --ner=20 --denoise --method=TV --print --random --weight
