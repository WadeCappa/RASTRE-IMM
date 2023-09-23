import json, sys

data = json.load(open(sys.argv[1], "r"))
timings = data["GranularRuntime_Milliseconds"]
world_size = data["WorldSize"]

diagnostics = [data["Levels"], data["BranchingFactor"]]
timing_output = [timings["Total"][0], timings["MPI_Gather"][0], timings["MPI_Gather"][world_size-1], timings["MaxKCover"][0], timings["MaxKCover"][world_size-1], timings["MPI_AllToAll"][0], timings["Sampling"][0]]
timing_output = list(map(lambda x: x / 1000, timing_output))

print(', '.join(list(map(lambda x: str(x), diagnostics + timing_output))))
