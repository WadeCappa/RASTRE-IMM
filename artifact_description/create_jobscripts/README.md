# Create jobscripts

The python scripts in this directory each generate a series of jobscripts. These scripts run individual jobs and are designed for use on the NERSC Perlmutter supercomputer. The job scripts generated, exhaustively match the experiments in Section 4 of the article. 

## Generating job scripts for different categories of experiments

For every table in section 4 the following python commands can be used to generate all the job scripts required to populate that table. These commands *do not nessessarly produce unique job scripts*, as there is overlap in the experiments across the tables.
* Table 4: Strong scaling performance of GreeDIMM for different inputs, varying $m$ up to $1024$ nodes. 
Command: `python3 strong_scaling_master.py <Account>`
* Table 5: Comparison of GreeDIMM with Ripples and IMM-mt. 
Command: `python3 comparison_master.py <Account>`
* Table 6: GreeDIMM-trunc: Runtime of the seed selection step (in sec) as a function of $\alpha$—the fraction of seeds sent from each sender—for Orkut with $m = 128$. 
Command: `python3 truncated_master.py <Account>`
* Table 7: Quality (measured in the expected number of activations
from the output seed set) as a function of $m$.
Command: `python3 quality_master.py <Account>`

These python scripts will output all generated experiment scripts into this directory.

## Submitting jobs on Perlmutter

One will be required to be part of a project on NERSC in order to submit jobs to the cluster. The command to submit a job is as follows: `sbatch <path-to-script>`.