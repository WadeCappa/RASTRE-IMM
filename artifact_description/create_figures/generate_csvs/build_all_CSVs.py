import csv_builder, sys

# need a function to output data in CSV format 
    # takes an array of valid keys as a header (determines order)
    # takes an array of tuples (ordered) containing a row label and a dictionary mapping keys to values
    # while building the talbe, only includes the vlaue in a row if it exists in the header.  

def main():
    builder = csv_builder.CSVBuilder()
    builder.build_strong_scaling("../../results/strong_scaling/", sys.argv[1] + "/strong_scaling.csv")

if __name__ == '__main__':
    main()