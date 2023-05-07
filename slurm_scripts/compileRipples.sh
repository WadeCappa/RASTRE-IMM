#!/bin/bash

module use /global/common/software/m3169/perlmutter/modulefiles

module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load openmpi

./waf configure --enable-mpi build_release
