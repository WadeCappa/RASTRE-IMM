import strong_scaling, imm_comparison, quality_drop_off, truncated_streaming
import sys

def main():
    strong_scalilng = strong_scaling.StrongScaling()
    strong_scalilng.build_strong_scaling("../../results/strong_scaling/", sys.argv[1] + "/strong_scaling.csv")

    comparison = imm_comparison.IMMComparison()
    comparison.build_comparison("../../results/imm/", "../../results/strong_scaling/", sys.argv[1])

    quality = quality_drop_off.QualityDropoff()
    quality.build_quality_csv("../../results/strong_scaling/", sys.argv[1] + "/streaming_quality_dropoff.csv")

    truncated = truncated_streaming.TruncatedStreaming()
    truncated.build_truncated_csv("../../results/truncated/orkut_small/", sys.argv[1] + "/truncated_results.csv")

if __name__ == '__main__':
    main()