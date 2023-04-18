import strong_scaling, imm_comparison, quality_drop_off
import sys

# need a function to output data in CSV format 
    # takes an array of valid keys as a header (determines order)
    # takes an array of tuples (ordered) containing a row label and a dictionary mapping keys to values
    # while building the talbe, only includes the vlaue in a row if it exists in the header.  

def main():
    strong_scalilng = strong_scaling.StrongScaling()
    strong_scalilng.build_strong_scaling("../../results/strong_scaling/", sys.argv[1] + "/strong_scaling.csv")

    comparison = imm_comparison.IMMComparison()
    comparison.build_comparison("../../results/imm/", "../../results/strong_scaling/", sys.argv[1])

    quality = quality_drop_off.QualityDropoff()
    quality.build_quality_csv("../../results/strong_scaling/", sys.argv[1] + "/streaming_quality_dropoff.csv")

if __name__ == '__main__':
    main()