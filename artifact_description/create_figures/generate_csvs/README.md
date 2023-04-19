# Generate CSVs

The scripts in this directory generate the tables in section 4 of the article as CSVs. These scripts rely on the raw experimental output files located in  `/artifact_description/results/`. 

**These python scripts must be run from this directory**. 

- `python3 build_all_CSVs.py <path/to/output>` will run all of the below python scripts and produce all of the CSVs described below.
- `python3 strong_scaling.py <path/to/output/>` will produce the CSV corresponding to table table 4 as `path/to/output/strong_scaling.csv`
- `python3 imm_comparison.py <path/to/output/>` will produce 4 tables corresponding to the data found in table 5 as 
  - `path/to/output/IC_comparison_quality.csv`
  - `path/to/output/IC_comparison_runtime.csv`
  - `path/to/output/LT_comparison_quality.csv`
  - `path/to/output/LT_comparison_runtime.csv`
- `python3 truncated_streaming.py <path/to/output/>` will produce the CSV corresponding to table 6 as `path/to/output/truncated_results.csv`
- `python3 quality_drop_off.py <path/to/output/>` will produce the CSV corresponding to table 7 `path/to/output/quality_dropoff.csv`
- `python3 scaling_breakdown.py <path/to/results/> <path/to/output>` will produce the 4 CSVs required by the `generate_plots` python scripts. These CSVs are
  - `path/to/output/wikipedia_timings.csv`
  - `path/to/output/wikipedia_receiver_timings.csv`
  - `path/to/output/livejournal_timings.csv`
  - `path/to/output/livejournal_receiver_timings.csv`