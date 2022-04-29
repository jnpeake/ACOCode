import os
import subprocess
import random
import csv
from unittest import result

ant_commands = ["./antmain.exe", "pcb442.tsp","500","16","40","99999", "20"]
num_seeds = 5
trials_per_seed = 1
num_tsp = 1
seeds = []
tsp = ["pcb442.tsp", "rat783.tsp", "pr1002.tsp", "pr2392.tsp"]
stagVals = [10,20,30,50,100,200,500,1000]

for i in range(num_seeds):
    seeds.append(str(random.randint(0,99999)))

with open('modifiednoveltytests16.csv', 'w') as csvfile:
    resultWriter = csv.writer(csvfile,
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for i in seeds:
        for j in stagVals:
            for k in tsp:
                resultWriter.writerow(["Seed: " + str(i) + " Stag: " + str(j) + " TSP: " + str(k)])
                ant_commands[5] = str(i)
                ant_commands[6] = str(j)
                ant_commands[1] = str(k)
                output = subprocess.check_output(ant_commands)
                    
                resultWriter.writerow(output.decode('utf-8').strip())