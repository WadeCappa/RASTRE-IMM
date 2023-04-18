#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -N 65
#SBATCH --ntasks-per-node=1
#SBATCH -J m65_IC_orkut_small
#SBATCH -o /Users/reetb/Documents/Research/GreedIMM/GreeDIMM/artifact_description/results/strong_scaling/orkut_small/m65_IC_orkut_small.o
#SBATCH -e /Users/reetb/Documents/Research/GreedIMM/GreeDIMM/artifact_description/results/strong_scaling/orkut_small/m65_IC_orkut_small.e


# module use /global/common/software/m3169/perlmutter/modulefiles
module use /global/common/software/m3169/perlmutter/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load cray-mpich
module load cray-libsci


srun -n 65 /Users/reetb/Documents/Research/GreedIMM/GreeDIMM/build/release/tools/mpi-greedi-im -i /Users/reetb/Documents/Research/GreedIMM/GreeDIMM/artifact_description/datasets/orkut_small_binary.txt -w -k 100 -p -u -d IC -e 0.13 -o /Users/reetb/Documents/Research/GreedIMM/GreeDIMM/artifact_description/results/strong_scaling/orkut_small/m65_IC_orkut_small.json --run-streaming=true --epsilon-2=0.077 --reload-binary