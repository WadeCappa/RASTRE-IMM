#!/bin/bash
module use /global/common/software/m3169/perlmutter/modulefiles



module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load cray-mpich
module load cray-libsci
# module load PrgEnv-nvidia
module load gpu/1.0
module load craype-accel-nvidia80
export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_ACCEL_TARGET=nvidia80
module load cudatoolkit

./waf configure --enable-mpi --enable-cuda build_release
