import csv_builder
import sys, math, json, os
from itertools import chain

class StrongScaling:
    def build_strong_scaling(self, parent_directory: str, output_file: str):
        builder = csv_builder.CSVBuilder()
        self.minimum_machines = sys.maxsize
        self.maximum_machines = 0

        def get_strong_scaling(parent_directory):
            experimental_map = dict() # maps datasets to a map of machines to runtime
            print(f"looking at {parent_directory}")

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
                    print(network, world_size, runtime)
                    self.minimum_machines = min(self.minimum_machines, world_size)
                    self.maximum_machines = max(self.maximum_machines, world_size)
                    if network not in experimental_map:
                        experimental_map[network] = dict()
                    experimental_map[network][world_size] = runtime

            return experimental_map

        experimental_map = get_strong_scaling(parent_directory)

        min_scale = int(math.log(self.minimum_machines - 1, 2))
        m_range = int(math.log(self.maximum_machines - 1, 2) - min_scale) + 1

        header = [int(math.pow(2, min_scale + x)) for x in range(min_scale + 1, m_range + 1)]
        print(header)
        builder.output_csv(output_file, header, experimental_map)

def main():
    strong_scalilng = StrongScaling()
    strong_scalilng.build_strong_scaling("../../results/strong_scaling/", sys.argv[1] + "/strong_scaling.csv")

if __name__ == '__main__':
    main()