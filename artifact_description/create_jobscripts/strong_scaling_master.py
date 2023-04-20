# Script to generate job scripts for all experiments in Table 4: Strong scaling performance of GreeDIMM for different inputs,
# varying m up to 1,024 nodes. 
# python3 strong_scaling_master.py <ACCOUNT>

import sys
import subprocess


nodes = [9, 17, 33, 65, 129, 257, 513, 1025]

datasets = ['Pokec', 'livejournal', 'orkut_small', 'orkut_big', 'wikipedia', 'friendster']

account = sys.argv[1]

for d in datasets:
	for n in nodes:
		if ((d == 'wikipedia' or d == 'friendster') and (n < 65)):
			continue
		else:
			if (d == 'orkut_small' or d == 'orkut_big' or d == 'wikipedia'):
				command = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 1'
			else:
				command = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 1'
			print(command)
			subprocess.call(command.split())