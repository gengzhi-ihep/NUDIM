#!/bin/bash
./projection_sim -i 4v6i.txt -p 128 --nx=128 --ny=128 --nz=128 --max=60 \
--min=-60 -t 3 -m 2D1D --noise --snr=1
