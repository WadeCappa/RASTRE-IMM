import csv_builder
import sys, os
from itertools import chain

class TruncatedStreaming:
    def __init__(self):
        pass

    def build_truncated_csv(self, single_network_directory: str, output: str):
        builder = csv_builder.CSVBuilder()
        headers = ['100', '80', '60', '40', '20']

        # get time
        output_files = builder.get_all_files(single_network_directory, lambda x: x.endswith(".json"))
        output_files = list(map(lambda x: {**{"Alpha": builder.get_alpha(x)}, **builder.get_data_from_json(x, ["Total", "WorldSize"])}, output_files))
        output_files = list(filter(lambda x: x['WorldSize'] == '129', output_files))
        output_rows = {"Time (in millisec)": {x["Alpha"]: x["Total"] for x in output_files}}

        # get quality
        quality_directories = builder.get_all_files(single_network_directory, lambda x: os.path.isdir(x))
        quality_files = list(chain.from_iterable([builder.get_all_files(quality_directory, lambda x: x.endswith(".json") and os.path.isfile(x)) for quality_directory in quality_directories]))
        quality_files = list(filter(lambda x: x["DiffusionModel"] == "IC", map(lambda json_file: {**{"Alpha": builder.get_alpha(json_file), "WorldSize": builder.get_machines(builder.strip_first_element(json_file))}, **builder.get_data_from_json(json_file, ["DiffusionModel", "Simulations"])}, quality_files)))
        quality_files = list(filter(lambda x: x['WorldSize'] == '129', quality_files))

        output_rows["Quality"] = dict()

        for q in quality_files:
            if q["Alpha"] in headers:
                output_rows["Quality"][q["Alpha"]] = builder.get_quality(q['Simulations'])

        first = None
        for header in headers:
            if first == None:
                first = output_rows["Quality"][header]
            output_rows["Quality"][header] = str((output_rows["Quality"][header] / first) * 100) + '%'

        output_rows = [(key, val) for key, val in output_rows.items()]
        for i, row in enumerate(output_rows):
            if row[0] == "Time (in millisec)":
                output_rows[i] = builder.convert_row_to_seconds(row)

        builder.output_csv(output, headers, output_rows)


def main():
    truncated_streaming = TruncatedStreaming()
    truncated_streaming.build_truncated_csv("../../results/truncated/orkut_small/", sys.argv[1] + "/truncated_results.csv")

if __name__ == '__main__':
    main()