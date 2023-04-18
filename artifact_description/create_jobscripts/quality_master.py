# Script to generate job scripts for all experiments in Table 7:Quality (measured in the expected number of activations from the output seed set) as a function of m.
# python3 quality_master.py <ACCOUNT>

import sys
import subprocess

nodes = [9, 17, 33, 65, 129]

datasets = ['github', 'HepPh', 'DBLP', 'Pokec', 'livejournal', 'orkut_small', 'orkut_big']

account = sys.argv[1]

for d in datasets:
	for n in nodes:
		if (d == 'github' or d == 'DBLP' or d == 'orkut_small' or d == 'orkut_big'):
			command = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 1'
		else:
			command = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 1'
		print(command)
		subprocess.call(command.split())