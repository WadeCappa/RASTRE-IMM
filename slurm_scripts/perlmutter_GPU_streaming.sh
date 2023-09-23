#!/bin/bash

#SBATCH -A m1641
#SBATCH -C gpu
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -N 5
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=4
#SBATCH -J m5_github_IC
#SBATCH -o /global/homes/r/reetb/cuda/results/jobs/github/m5_github_IC.o
#SBATCH -e /global/homes/r/reetb/cuda/results/jobs/github/m5_github_IC.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wade.cappa@wsu.edu

module use /global/common/software/m3169/perlmutter/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# module load PrgEnv-nvidia
module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load cray-mpich
module load cray-libsci
module load gpu/1.0
module load craype-accel-nvidia80
export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_ACCEL_TARGET=nvidia80
module load cudatoolkit




srun -n 5 ./build/release/tools/mpi-greedimm -i /global/cfs/cdirs/m1641/network-data/Binaries/github_IC_binary.txt --streaming-gpu-workers 4 -w -k 100 -p -d IC -e 0.13 -o /global/homes/r/reetb/cuda/results/jobs/github/m5_github_IC.json --run-streaming=true --epsilon-2=0.077 --reload-binary -u
