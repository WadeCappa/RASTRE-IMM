import sys, os

runtime = dict() # maps datasets to a map of machines to runtime

for file in os.listdir(sys.argv[1]):
    if os.path.isfile(file):
        continue
    for network_directory in os.listdir(file):