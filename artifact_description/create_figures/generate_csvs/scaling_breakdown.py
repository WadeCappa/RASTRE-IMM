import csv_builder
import sys, os, math
from itertools import chain

        # --- SHARED --- 
        # Samping time: 35367.6
        # AlltoAll time: 37826.1
        # Receive Broadcast: 205.391
        # --- SENDER --- 
        # Select Next Seed: 1489.78
        # Send Next Seed: 81.3104
        # Total Send Time: 3220.8
        # --- RECEIVER --- 
        # Initialize Buckets: 0
        # Receive Next Seed: 0
        # Insert Into Buckets: 0
        # Handling received data (inserting into matrix and copying from buffer): 0
        # Atomic Update (receiver side): 0
        # Total Global Streaming Time: 0
        # exited
        # [2023-04-01 20:47:28.025] [console] [info] IMM MPI+OpenMP+CUDA : 77368.377733ms
        # [2023-04-01 20:47:28.026] [console] [info] IMM World Size : 9
        # [2023-04-01 20:47:28.026] [console] [info] IMM Rank : 6

class ScalingBreakdown:
    def __init__(self):
        self.DELIMITER = ' --- SHARED --- \n'
        self.RECEIVER = 'receiver'
        self.SENDER = 'sender'
        self.SHARED = 'shared'
        self.SAMP = 'Samping time'
        self.ALL2ALL = 'AlltoAll time'
        self.BCAST = 'Receive Broadcast'
        self.SELECT_SEED = 'Select Next Seed'
        self.SEND_SEED = 'Send Next Seed'
        self.TOTAL_SEND = 'Total Send Time'
        self.INIT_BUCKET = 'Initialize Buckets'
        self.RECV_SEED = 'Receive Next Seed'
        self.INSERT = 'Insert Into Buckets'
        self.HANDLE = 'Handling received data'
        self.ATOMIC = 'Atomic Update'
        self.TOTAL_RECV = 'Total Global Streaming Time'
        self.TOTAL = 'IMM MPI+OpenMP+CUDA'

        self.KEYWORDS = [
            self.SAMP, 
            self.ALL2ALL, 
            self.BCAST, 
            self.SELECT_SEED, 
            self.SEND_SEED,
            self.TOTAL_SEND, 
            self.INIT_BUCKET, 
            self.RECV_SEED, 
            self.INSERT, 
            self.HANDLE, 
            self.ATOMIC, 
            self.TOTAL_RECV,
            self.TOTAL
        ]

        self.SENDER_KEYWORDS = set([            
            self.SELECT_SEED, 
            self.SEND_SEED,
            self.TOTAL_SEND,
            self.SAMP,
            self.ALL2ALL,
            self.BCAST
        ])

        self.RECEIVER_KEYWORDS = set([
            self.INIT_BUCKET, 
            self.RECV_SEED, 
            self.INSERT, 
            self.HANDLE, 
            self.ATOMIC, 
            self.TOTAL_RECV,
            self.SAMP,
            self.ALL2ALL,
            self.BCAST
        ])

        self.SHARED_KEYWORDS= set({
            self.TOTAL
        })


    def load_output_file(self, output_file: str):
        with open(output_file, 'r') as stream:
            return stream.read()
        
    def get_timings_from_block(self, block: list[str]):
        res = {self.RECEIVER: {key: 0 for key in self.RECEIVER_KEYWORDS}, self.SENDER: {key: 0 for key in self.SENDER_KEYWORDS}, self.SHARED: {key: 0 for key in self.SHARED_KEYWORDS}}
        for line in block:
            for key in self.KEYWORDS:
                if key in line:
                    timing = float(line.split(': ')[-1]) if key != self.TOTAL else float(line.split(': ')[-1][:-2])
                    # res[self.RECEIVER if key in self.RECEIVER_KEYWORDS else (self.SENDER if key in self.SENDER_KEYWORDS else self.SHARED)][key] = timing
                    # print(key, key in self.SHARED_KEYWORDS, key in self.RECEIVER_KEYWORDS, key in self.SENDER_KEYWORDS)
                    if key in self.SHARED_KEYWORDS:
                        res[self.SHARED][key] = timing 
                    if key in self.RECEIVER_KEYWORDS:
                        res[self.RECEIVER][key] = timing
                    if key in self.SENDER_KEYWORDS:
                        res[self.SENDER][key] = timing

    
        return res
    
    def block_is_receiver(self, block: dict[str, float]):
        return self.TOTAL_RECV in block[self.RECEIVER] and block[self.RECEIVER][self.TOTAL_RECV] != 0
    
    def block_is_sender(self, block: dict[str, float]):
        return not self.block_is_receiver(block) and self.TOTAL_SEND in block[self.SENDER]

    def deconstruct_ouput_file(self, output_data: str) -> dict[str, dict[str, str]]: # returns receiver or sender -> dict of time_lable -> time
        file_data = {self.RECEIVER: {key: 0 for key in self.RECEIVER_KEYWORDS}, self.SENDER: {key: 0 for key in self.SENDER_KEYWORDS}, self.SHARED: {key: 0 for key in self.SHARED_KEYWORDS}}

        blocks = [self.get_timings_from_block(block.split('\n')) for block in output_data.split(self.DELIMITER)]\
        
        best_time = 0

        for block in blocks:
            if self.block_is_receiver(block):
                for key in self.RECEIVER_KEYWORDS:
                    file_data[self.RECEIVER][key] = block[self.RECEIVER][key]
                file_data[self.SHARED][self.TOTAL] = block[self.SHARED][self.TOTAL]
                
            elif self.block_is_sender(block):
                if block[self.SHARED][self.TOTAL] > best_time:
                    best_time = block[self.SHARED][self.TOTAL]
                    for key in self.SENDER_KEYWORDS:
                        file_data[self.SENDER][key] = block[self.SENDER][key]

        return file_data
        

    def build_breakdown_csv(self, single_network_directory: str, machine_range: list[str], output: str):
        builder = csv_builder.CSVBuilder()
        headers = ["Total", "S-Sampling", "S-All2All", "S-SeedSelect", "R-Sampling", "R-All2All", "R-SeedSelect"]

        output_files = builder.get_all_files(single_network_directory, lambda x: x.endswith(".o"))
        output_files = list(filter(lambda x: "LT" not in x.split("/")[-1].split('_') and str(int(builder.get_machines(x)) - 1) in machine_range, output_files))

        trials = [{**self.deconstruct_ouput_file(self.load_output_file(output_file)), **{'WorldSize': str(int(builder.get_machines(output_file))-1)}} for output_file in output_files]
        # print(trials)

        rows = [builder.convert_row_to_seconds((trial['WorldSize'], {
            headers[0]: trial[self.SHARED][self.TOTAL], 
            headers[1]: trial[self.SENDER][self.SAMP], 
            headers[2]: trial[self.SENDER][self.ALL2ALL], 
            headers[3]: trial[self.SENDER][self.TOTAL_SEND], 
            headers[4]: trial[self.RECEIVER][self.SAMP],
            headers[5]: trial[self.RECEIVER][self.ALL2ALL], 
            headers[6]: trial[self.RECEIVER][self.TOTAL_RECV], 
        })) for trial in trials]

        builder.output_csv(output, headers, rows)

    def build_receiver_breakdown_csv(self, single_network_directory: str, machine_range: list[str], output: str):
        builder = csv_builder.CSVBuilder()
        headers = ["Total", "Waiting_and_Processing", "Communication", "Push_to_buckets"]

        output_files = builder.get_all_files(single_network_directory, lambda x: x.endswith(".o"))
        output_files = list(filter(lambda x: "LT" not in x.split("/")[-1].split('_') and str(int(builder.get_machines(x)) - 1) in machine_range, output_files))

        trials = [{**self.deconstruct_ouput_file(self.load_output_file(output_file)), **{'WorldSize': str(int(builder.get_machines(output_file))-1)}} for output_file in output_files]
        # print(trials)

        rows = [builder.convert_row_to_seconds((trial['WorldSize'], {
            headers[0]: trial[self.RECEIVER][self.TOTAL_RECV], 
            headers[1]: trial[self.RECEIVER][self.INSERT], 
            headers[2]: trial[self.RECEIVER][self.RECV_SEED], 
            headers[3]: trial[self.RECEIVER][self.HANDLE]
        })) for trial in trials]

        builder.output_csv(output, headers, rows)


def main():
    scaling_breakdown = ScalingBreakdown()
    scaling_breakdown.build_breakdown_csv(sys.argv[1] + "/wikipedia", [str(int(math.pow(2, x))) for x in range(6, 10 + 1)], sys.argv[2] + "/wikipedia_timings.csv")
    scaling_breakdown.build_breakdown_csv(sys.argv[1] + "/livejournal", [str(int(math.pow(2, x))) for x in range(4, 10 + 1)], sys.argv[2] + "/livejournal_timings.csv")

    scaling_breakdown.build_receiver_breakdown_csv(sys.argv[1] + "/wikipedia", [str(int(math.pow(2, x))) for x in range(6, 10 + 1)], sys.argv[2] + "/wikipedia_receiver_timings.csv")
    scaling_breakdown.build_receiver_breakdown_csv(sys.argv[1] + "/livejournal", [str(int(math.pow(2, x))) for x in range(4, 10 + 1)], sys.argv[2] + "/livejournal_receiver_timings.csv")

if __name__ == '__main__':
    main()