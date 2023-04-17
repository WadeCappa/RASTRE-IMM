import sys, os, json, math
from itertools import chain

class CSVBuilder:
    def __init__(self):
        pass

    def build_comparison(self, directory, output):
        # two runs, one for "IC" one for "LT"
            # open directory, scan through all files, get subdirectories
        pass

    def get_all_files(self, directory, filter_func):
        return list(filter(filter_func, list(map(lambda x: directory + "/" + x, os.listdir(directory)))))

    def build_strong_scaling(self, parent_directory, output_file):
        self.minimum_machines = sys.maxsize
        self.maximum_machines = 0

        def get_strong_scaling(parent_directory):
            experimental_map = dict() # maps datasets to a map of machines to runtime
            print(f"looking at {parent_directory}")

            key_directories = self.get_all_files(parent_directory, lambda x: os.path.isdir(x))
            outputs = list(chain.from_iterable([self.get_all_files(directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for directory in key_directories]))

            for result in outputs:
                input_stream = open(result)
                data = json.load(input_stream)
                input_stream.close()

                network =  data[0]["Input"].split("/")[-1]
                world_size = data[0]["WorldSize"] if "WorldSize" in data[0] and data[0]["WorldSize"] % 2 == 1 else None
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
        num_networks = len(experimental_map)
        print(m_range, num_networks, self.maximum_machines, self.minimum_machines, int(math.log(self.maximum_machines - 1, 2)), min_scale)

        with open(output_file, "w") as out:
            out.write(", " + ', '.join([str(int(math.pow(2, min_scale + x))) for x in range(m_range)]) + '\n')
            for network, results in experimental_map.items():
                row = [None for _ in range(int(m_range))]
                for world_size, time in results.items():
                    row[int(math.log(world_size-1, 2) - min_scale)] = time

                out.write(network.replace("_binary.txt", "") + ', ' + ', '.join(list(map(lambda x: str(x) if x != None else "", row))) + "\n")

if __name__ == '__main__':
    builder = CSVBuilder()
    builder.build_strong_scaling(sys.argv[1], sys.argv[2])
