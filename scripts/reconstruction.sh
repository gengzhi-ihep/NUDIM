#!/bin/bash
./reconstruction --prefix=prj --pad=128 --nx=128 --nz=128 --max=60 \
--min=-60 -t 3 --padangle=1 --nhio=0 --ner=20 --denoise --method=TV --print --random --weight
