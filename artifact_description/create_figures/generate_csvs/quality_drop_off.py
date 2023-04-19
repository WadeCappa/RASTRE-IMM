import csv_builder
import sys, os, math
from itertools import chain

class QualityDropoff:
    def __init__(self):
        pass 

    def build_quality_csv(self, network_directory: str, output_file: str):
        builder = csv_builder.CSVBuilder()
        headers = [int(math.pow(2, x)) for x in range(int(math.log(8, 2)), int(math.log(128, 2) + 1))]
        networks = ["github", "HepPh", "orkut_small", "orkut_big", "DBLP", "Pokec", "livejournal"]
        
        network_directories = builder.get_all_files(network_directory, lambda x: os.path.isdir(x))
        quality_directories = list(chain.from_iterable([builder.get_all_files(possible_quality_directory, lambda x: os.path.isdir(x)) for possible_quality_directory in network_directories]))
        quality_data = list(chain.from_iterable([builder.get_all_files(quality_directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for quality_directory in quality_directories]))
        quality_data = list(filter(lambda x: x[1] != None and x[1]["DiffusionModel"] == "IC", map(lambda json_file: [json_file, builder.get_data_from_json(json_file, ["DiffusionModel", "Input", "Simulations"])], quality_data)))
        quality_data = list(filter(lambda x: x['WorldSize'] in headers, [{**{"WorldSize": int(data[0].split("/")[-1].split("_")[0].replace('m', '')) - 1 if (data[0].split("/")[-1].split("_")[0].replace('m', '')).isdigit() else 0}, **data[1]} for data in quality_data]))

        rows = {key: dict() for key in networks} # network -> dict(WorldSize -> quality_percentage) where WorldSize in headers

        for experiment in quality_data:
            network = builder.get_network(experiment['Input'])
            if network == "githubSmall":
                network = "github"
            if network in rows and experiment['WorldSize'] in headers:
                rows[network][experiment['WorldSize']] = builder.get_quality(experiment["Simulations"])

        for key, row in rows.items():
            first = None
            for header in headers:
                if first == None:
                    first = row[header]
                row[header] = str((row[header] / first) * 100) + '%'
                
        
        builder.output_csv(output_file, headers, [(key, val) for key, val in rows.items()])

def main():
    quality_dropoff = QualityDropoff()
    quality_dropoff.build_quality_csv("../../results/strong_scaling/", sys.argv[1] + "/quality_dropoff.csv")

if __name__ == '__main__':
    main()