import csv_builder
import sys, math, json, os
from itertools import chain

class StrongScaling:
    def build_strong_scaling(self, parent_directory: str, output_file: str):
        builder = csv_builder.CSVBuilder()
        self.minimum_machines = sys.maxsize
        self.maximum_machines = 0

        networks = ["Pokec", "livejournal", "orkut_small", "orkut_big", "wikipedia", "friendster"]

        def get_strong_scaling(parent_directory):
            experimental_map = {key: dict() for key in networks} # maps datasets to a map of machines to runtime

            key_directories = builder.get_all_files(parent_directory, lambda x: os.path.isdir(x))
            outputs = list(chain.from_iterable([builder.get_all_files(directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for directory in key_directories]))

            for result in outputs:
                input_stream = open(result)
                data = json.load(input_stream)
                input_stream.close()

                network = data[0]["Input"].split("/")[-1]
                world_size = data[0]["WorldSize"] - 1 if "WorldSize" in data[0] and data[0]["WorldSize"] % 2 == 1 else None
                runtime = data[0]["Total"] if "Total" in data[0] else None

                if world_size != None and runtime != None and network != None and "DiffusionModel" in data[0] and data[0]["DiffusionModel"] == "IC":
                    self.minimum_machines = min(self.minimum_machines, world_size)
                    self.maximum_machines = max(self.maximum_machines, world_size)
                    network_name = builder.get_network(data[0]["Input"])
                    if network_name in experimental_map:
                        experimental_map[network_name][world_size] = runtime

            return experimental_map

        experimental_map = get_strong_scaling(parent_directory)

        min_scale = int(math.log(self.minimum_machines - 1, 2))
        # m_range = int(math.log(self.maximum_machines - 1, 2) - min_scale) + 1

        header = [int(math.pow(2, min_scale + x)) for x in range(int(math.log(8, 2)), int(math.log(1024, 2)) + 1)]
        builder.output_csv(output_file, header, [builder.convert_row_to_seconds((key, val)) for key, val in experimental_map.items()])

def main():
    strong_scalilng = StrongScaling()
    strong_scalilng.build_strong_scaling("../../results/strong_scaling/", sys.argv[1] + "/strong_scaling.csv")

if __name__ == '__main__':
    main()