#!/bin/bash
/usr/lib64/openmpi/bin/mpirun -n 9 ../../bin/reconstruction_parallel --prefix=prj --pad=512 --nx=512 --nz=512 \
--max=63 --min=-55 -t 2 --padangle=1 --nhio=0 --ner=20 --denoise --method=TV --print --random --weight
