# Datasets

The scripts in this directory can be used to prepare the input datasets from the SNAP and Konect database. 

* Download the datasets (only the lines corresponding to the dataset required can be left uncommented) into this directory.
Command:
`./get_datasets.sh`

* Compile the datasets (only the lines corresponding to the dataset required can be left uncommented) to generate the binaries. This will dump the datasets in the format format after generating edge weights uniformaly distributed between $(0,0.1)$.
Command:
`./compile_datasets.sh`
