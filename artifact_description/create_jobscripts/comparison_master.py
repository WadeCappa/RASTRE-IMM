# Script to generate job scripts for all experiments in Table 5: Comparison of GreeDIMM with Ripples and IMM-mt. 
# python3 comparison_master.py <ACCOUNT>

import sys
import subprocess

nodes = [1, 513]

datasets = ['github', 'HepPh', 'DBLP', 'Pokec', 'livejournal', 'orkut_small', 'orkut_big', 'wikipedia', 'friendster']

account = sys.argv[1]

for d in datasets:
	for n in nodes:
		if (n == 513):
			if (d == 'github' or d == 'DBLP' or d == 'orkut_small' or d == 'orkut_big' or d == 'wikipedia'):
				command_IC = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 1'
				command_LT = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 0'
				command_IC_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n - 1) + ' ' + str(d) + ' 1' + ' 1'
				command_LT_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n - 1) + ' ' + str(d) + ' 1' + ' 0'
			else:
				command_IC = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 1'
				command_LT = 'python3 strong_scaling_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 0'
				command_IC_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n - 1) + ' ' + str(d) + ' 0' + ' 1'
				command_LT_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n - 1) + ' ' + str(d) + ' 0' + ' 0'
			
			subprocess.call(command_IC.split())
			subprocess.call(command_LT.split())
			subprocess.call(command_IC_IMM.split())
			subprocess.call(command_LT_IMM.split())
		else:
			if (d == 'github' or d == 'DBLP' or d == 'orkut_small' or d == 'orkut_big' or d == 'wikipedia'):
				command_IC_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 1'
				command_LT_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 1' + ' 0'
			else:
				command_IC_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 1'
				command_LT_IMM = 'python3 imm_jobs.py ' + str(account) + ' ' + str(n) + ' ' + str(d) + ' 0' + ' 0'
		
			subprocess.call(command_IC_IMM.split())
			subprocess.call(command_LT_IMM.split())