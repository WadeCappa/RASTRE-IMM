import os, functools, json

class CSVBuilder:
    def __init__(self):
        pass

    def get_network(self, input_path):
        return input_path.split("/")[-1].replace("_binary.txt", '')
    
    def get_alpha(self, file_path):
        return file_path.split('/')[-1].split('_')[0].replace('a', '')
    
    def get_machines(self, file_path):
        return file_path.split("/")[-1].split("_")[0].replace('m', '')
    
    def strip_first_element(self, file_path):
        path = file_path.split('/')
        path[-1] = '_'.join(path[-1].split('_')[1:])
        return '/'.join(path)
    
    def get_quality(self, simulations):
        simulations = json.loads('{ "data":' + simulations + '}')['data']
        # print(simulations)
        return functools.reduce(lambda x, y: x + y, list(map(lambda x: int(x[0]), simulations))) / len(simulations)

    def load_json(self, json_file: str) -> dict:
        if os.path.getsize(json_file) == 0:
            return dict()
        input_stream = open(json_file)
        json_data = json.load(input_stream)
        input_stream.close()
        return json_data[0]

    def does_json_conform(self, json_data: dict, data: list[str]) -> bool:
        return functools.reduce(lambda x, y: x and y, [val in json_data for val in data])

    def get_data_from_json(self, json_file: str, data: list[str]) -> dict:
        json_data = self.load_json(json_file)
        if not self.does_json_conform(json_data, data):
            return None 
        res = dict()
        for val in data:
            res[val] = str(json_data[val])
        return res

    def get_all_files(self, directory, filter_func):
        return list(filter(filter_func, list(map(lambda x: directory + "/" + x, os.listdir(directory)))))
    
    def convert_row_to_seconds(self, row: tuple[str, dict[str, float]]):
        return (row[0], {key: float(val) / 1000 for key, val in row[1].items()})
    
    def output_csv(self, outfile: str, header: list, rows: list[list[str, dict]]):
        with open(outfile, "w") as out:
            out.write('\t' + '\t'.join(list(map(lambda x: str(x), header))) + '\n')
            for network, row in rows:
                new_row = []       

                for val in header:
                    new_row.append(str(row[val]) if val in row else "")  
                      
                out.write(network + '\t' + '\t'.join(new_row) + '\n')