import csv_builder
import sys, os
from itertools import chain

class IMMComparison:
    def __init__(self):
        pass

    def build_comparison(self, imm_directory: str, greedimm_directory, output: str):
        # two runs, one for "IC" one for "LT"
            # open imm_directory, scan through all files, get subdirectories
                # for every subdirectory, record the diffusion_model_1 and diffusion_model_256 run
                # record the diffusion_model_1 and diffusion_model_256 quality
            # open greedimm_directory, scan through all files get subdirectories
                # for every subdirectory record difusion_model_257 run 
                # record the diffusion_model 257 quality
            # build header, output as table_diffusion_model

        builder = csv_builder.CSVBuilder()

        for diffusion_model in ["IC", "LT"]:
            imm_network_directories = builder.get_all_files(imm_directory, lambda x: os.path.isdir(x))
            imm_output_times = list(chain.from_iterable([builder.get_all_files(experimental_result, lambda x: x.endswith(".json") and os.path.isfile(x)) for experimental_result in imm_network_directories]))
            imm_output_times = list(map(lambda file: builder.get_data_from_json(file, ["DiffusionModel", "Input", "WorldSize", "Total"]), imm_output_times))
            imm_output_times = list(filter(lambda x: x != None and (x['WorldSize'] == '1' or x['WorldSize'] == '256') and x['DiffusionModel'] == diffusion_model, imm_output_times))
            # print(imm_output_times)

            imm_quality_directories = list(chain.from_iterable([builder.get_all_files(possible_quality_directory, lambda x: os.path.isdir(x)) for possible_quality_directory in imm_network_directories]))
            imm_quality_data = list(chain.from_iterable([builder.get_all_files(quality_directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for quality_directory in imm_quality_directories]))
            imm_quality_data = list(filter(lambda x: x[1] != None and x[1]["DiffusionModel"] == diffusion_model, map(lambda json_file: [json_file, builder.get_data_from_json(json_file, ["DiffusionModel", "Input", "Simulations"])], imm_quality_data)))
            imm_quality_data = list(filter(lambda x: x['WorldSize'] == '1' or x['WorldSize'] == '256', [{**{"WorldSize": data[0].split("/")[-1].split("_")[1]}, **data[1]} for data in imm_quality_data]))
            # print(imm_quality_data)

            greedimm_network_directories = builder.get_all_files(greedimm_directory, lambda x: os.path.isdir(x))
            greedimm_output_times = list(chain.from_iterable([builder.get_all_files(experimental_result, lambda x: x.endswith(".json") and os.path.isfile(x)) for experimental_result in greedimm_network_directories]))
            greedimm_output_times = list(map(lambda file: builder.get_data_from_json(file, ["DiffusionModel", "Input", "WorldSize", "Total"]), greedimm_output_times))
            greedimm_output_times = list(filter(lambda x: x != None and x['WorldSize'] == '257' and x['DiffusionModel'] == diffusion_model, greedimm_output_times))

            greedimm_quality_directories = list(chain.from_iterable([builder.get_all_files(possible_quality_directory, lambda x: os.path.isdir(x)) for possible_quality_directory in greedimm_network_directories]))
            greedimm_quality_data = list(chain.from_iterable([builder.get_all_files(quality_directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for quality_directory in greedimm_quality_directories]))
            greedimm_quality_data = list(filter(lambda x: x[1] != None and x[1]["DiffusionModel"] == diffusion_model, map(lambda json_file: [json_file, builder.get_data_from_json(json_file, ["DiffusionModel", "Input", "Simulations"])], greedimm_quality_data)))
            greedimm_quality_data = list(filter(lambda x: x['WorldSize'] == '257', [{**{"WorldSize": data[0].split("/")[-1].split("_")[0].replace('m', '')}, **data[1]} for data in greedimm_quality_data]))

            # print(greedimm_quality_data)
            # print(greedimm_output_times)

            headers = ["IMM-mt", "Ripples", "GreeDIMM"]

            def get_header(data: dict) -> str:
                machines = data['WorldSize']
                if machines == '1': return headers[0]
                if machines == '256': return headers[1]
                if machines == '257': return headers[2]

            quality_rows = dict() # networkname -> (header -> quality)
            runtime_rows = dict() # networkname -> (header -> runtime)

            for data_set in [greedimm_output_times, imm_output_times]:
                for data in data_set:
                    network = builder.get_network(data["Input"])
                    if network not in runtime_rows:
                        runtime_rows[network] = dict()
                    runtime_rows[network][get_header(data)] = data["Total"]

            for data_set in [greedimm_quality_data, imm_quality_data]:
                for data in data_set:
                    network = builder.get_network(data["Input"])
                    if network not in quality_rows:
                        quality_rows[network] = dict()
                    quality_rows[network][get_header(data)] = builder.get_quality(data["Simulations"])

            # print(quality_rows)
            builder.output_csv(output + f"/{diffusion_model}_comparison_runtime.csv", headers, runtime_rows)
            builder.output_csv(output + f"/{diffusion_model}_comparison_quality.csv", headers, quality_rows)


            # rows will be a dict of network names to a nested dictionary
                # this nested dictionary will map an algoirthm name (same as in header) to runtime or quality
            # headers will be a list of each algorithm name
            # 4 tables will be produced, LT for runtime and quality, IC for runtime and quality

            # now that you have all the correct data, build 4 tables.
                # to do this, first add a new string to all dictionaries of the algorithm used. If machines == 1, IMM-mt, else Ripples for the IMM jobs. 
                # The headers will just match the types of algorithms used. Combine IMM and ripples datasets

def main():
    imm_comparison = IMMComparison()
    imm_comparison.build_comparison("../../results/imm/", "../../results/strong_scaling/", sys.argv[1])

if __name__ == '__main__':
    main()