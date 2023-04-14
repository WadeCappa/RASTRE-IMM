import os
import sys
import json

for filename in os.listdir(sys.argv[1]):
    data = open(sys.argv[1] + "/" + filename)
    q = json.load(data) 
    data.close()

    result = sum(list(map(lambda x: x[0], q[0]["Simulations"]))) / 5
    print(f"{result} from {filename}")
