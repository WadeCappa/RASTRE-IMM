import sys
import os


account = sys.argv[1]
nodes = sys.argv[2]
dataset = sys.argv[3]
directed = '-u' if (sys.argv[4] == '1') else ''
model = 'IC' if (sys.argv[5] == '1') else 'LT'
output_path = os.path.abspath('../results/imm/' + dataset + '/imm_' + nodes + '_' + model + '_' + dataset)
tool_path = os.path.abspath('../../build/release/tools/mpi-imm')
dataset_path = os.path.abspath('../datasets/' + dataset + '_binary.txt')


script = '''#!/bin/bash

#SBATCH -A {}
#SBATCH -C cpu
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -N {}
#SBATCH --ntasks-per-node=1
#SBATCH -J imm_{}_{}_{}
#SBATCH -o {}.o
#SBATCH -e {}.e


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


srun -n {} {} -i {} -w -k 100 -p {} -d {} -e 0.13 -o {}.json --reload-binary'''.format(account, nodes, nodes, model, dataset, output_path, output_path, nodes, tool_path, dataset_path, directed, model, output_path)


with open('imm_' + str(nodes) + '_' + model + '_' + dataset + '.sh', 'w') as f:
    f.write(script)