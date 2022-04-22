# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 09:02:06 2022

@author: 55125529
"""

import pandas as pd

class Result():
    def __init__(self, seed, tsp, stag, time, qual):
        self.seed = seed
        self.tsp = tsp
        self.stag = stag
        self.time = time
        self.qual = qual
        
    def to_dict(self):
        return {
            "seed": self.seed,
            "tsp": self.tsp,
            "stag": int(self.stag),
            "time": float(self.time),
            "qual": int(float(self.qual)),
            }
        
f = open("noveltytests8.csv", "r")
results = []
count = 0
for index, line in enumerate(f):
    if line[0] == "S":
        splitLine = line.split()
        currSeed = splitLine[1]
        currStag = splitLine[3]
        currTSP = splitLine[5]
    if line[0] == "s":
        splitLine = line.split()
        currQual = splitLine[2]
        currTime = splitLine[4]
    count += 1

    if count == 4:
        results.append(Result(currSeed,currTSP,currStag,currTime,currQual))
        count = 0

df = pd.DataFrame.from_records([r.to_dict() for r in results])

df = df.sort_values(by = ["tsp", "stag"])
pcbDF = df.loc[df['tsp'] == "pcb442.tsp"]
ratDF = df.loc[df['tsp'] == "rat783.tsp"]
prDF = df.loc[df['tsp'] == "pr1002.tsp"]

pcbMean = (pcbDF.groupby(['stag']).mean()).reset_index()
ratMean = (ratDF.groupby(['stag']).mean()).reset_index()
prMean = (prDF.groupby(['stag']).mean()).reset_index()

print(pcbMean.columns)
pcbMean.plot(x = 'stag', y = 'qual', marker = "o")
ratMean.plot(x = 'stag', y = 'qual', marker = "o")
prMean.plot(x = 'stag', y = 'qual', marker = "o")