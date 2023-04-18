# Create jobscripts

The scripts in this directory have been orchestrated to generate the job scripts required to run all the experiments in Section 4 of the article on the NERSC Perlmutter supercomputer. 


### Submitting jobs on Perlmutter

One will require to be part of a project on NERSC Jobs in order to submit jobs to the cluster. The command to submit a job is as follows: 
<!-- TODO: fix which directory to call from (Is this even required since paths are absolute?)-->
`sbatch <path-to-script>`

### Scripts for different categories of experiments
* Table 4: Strong scaling performance of GreeDIMM for different inputs, varying $m$ up to $1024$ nodes. 
Command: 
`python3 strong_scaling_master.py <Account>`
* Table 5: Comparison of GreeDIMM with Ripples and IMM-mt. 
Command:
`python3 comparison_master.py <Account>`
* Table 6: GreeDIMM-trunc: Runtime of the seed selection step (in sec) as a function of $\alpha$—the fraction of seeds sent from each sender—for Orkut with $m = 128$. 
Command:
`python3 truncated_master.py <Account>`
* Table 7: Quality (measured in the expected number of activations
from the output seed set) as a function of $m$.
Command:
`python3 quality_master.py <Account>`






