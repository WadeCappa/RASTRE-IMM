# Datasets

The scripts in this directory can be used to prepare the input datasets used throughout the article. These datasets come from the SNAP and Konect databases. 

* The `get_datasets.sh` bash script will download all 9 datasets used in the article. If all of the datasets are not required each dataset can be individually downloaded using the corresponding commands in the script. 

* The `compile_datasets.sh` script produces all the inputs used by the `GreeDIMM` application. These inputs are generated such that each edge is assigned a weight uniformly at random in the range $(0,0.1)$.