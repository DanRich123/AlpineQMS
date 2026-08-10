#!/bin/bash

# optional flags
OPT_FLAGS="-O3 -qmkl"

# load modules
# ml purge # useful if many packages have been loaded and we want to start fresh
ml intel
ml mkl

# Now compile my code
ifx ${OPT_FLAGS} QMS_3D.f90 -o QMS_3D