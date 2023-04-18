# Script to generate job scripts for all experiments in Table 6: GreeDIMM-trunc: Runtime of the seed selection step (in sec)
# as a function of 𝛼—the fraction of seeds sent from each sender for Orkut with m = 128.
# python3 truncated_master.py <ACCOUNT>

import sys
import subprocess


nodes = [129]

datasets = ['orkut_small']

alpha = [100, 80, 60, 40 , 20]

account = sys.argv[1]

for d in datasets:
	for n in nodes:
		for a in alpha:
			command = 'python3 truncated_jobs.py ' +  str(account) + ' ' + str(n) + ' ' +  str(d) + ' 0 1 ' + str(a)
			print(command)
			subprocess.call(command.split())