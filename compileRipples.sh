#!/bin/bash

module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.22.1
module load openmpi/4.1.2

./waf configure --enable-mpi build_release
