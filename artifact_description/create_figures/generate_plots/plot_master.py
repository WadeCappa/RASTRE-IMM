
import sys
import subprocess

datasets = ['livejournal', 'wikipedia']

for d in datasets:
	command1 = 'python3 greedi_plots.py --data ' + d + '_data.csv --output Figure4_' + d + '_breakdown.PDF' 
	command2 = 'python3 greedi_plots_receiver.py --data ' + d + '_receiver_data.csv --output Figure5_' + d + '_receiver_breakdown.PDF' 
	subprocess.call(command1.split())
	subprocess.call(command2.split())