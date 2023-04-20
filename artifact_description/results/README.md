# Results

This directory contains the experimental outputs generated from all the runs conducted on the NERSC Perlmutter supercomputer. Scripts to parse the logs to generate the CSVs and PDFs corresponding to the tables and plots respectively in Section 4 can be found in the directory `/artifact_description/create_figures/`.

The outputs for each experiment are split into three distinct files
- The `<experiment>.json` file contains the output of an experiment. 
- The `<experiment>.o` file contains the diagnostic information collected during the runtime of the experiment. 
- The `<experiment>.e` is populated if an error occures during the runtime of the experiment. 

The experimental results used in the article are divided into three subdirectories;
- `/imm/` contains the experimental results for the `IMM-mt` and `Ripples` columns in table 5. 
- `/strong_scaling/` contains the experimental results for 
  - table 4
  - the results for the `GreeDIMM` column in table 5
  - table 7
  - figure 4 and figure 5
- `/truncated/` contains the experimental results for table 6

## Quality
To recalculate the quality estimation for a completed job the `./get_quality.sh <path/to/GreeDIMM> <path/to/experiment/network/directory/> <path/to/input/dataset/>` script can be used. This iterates through all of the `*.json` files in the `<path/to/experiment/network/directory/>` directory and outputs the quality data into the `<path/to/experiment/network/directory/quality/>` directory.