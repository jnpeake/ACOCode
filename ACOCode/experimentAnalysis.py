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
    
def readFile(inputFile):
    count = 0
    tempResults = []
    for index, line in enumerate(inputFile):
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

        if count == 2:
            tempResults.append(Result(currSeed,currTSP,currStag,currTime,currQual))
            count = 0
    return tempResults
    
        
faco = open("acotests16.csv", "r")
fnovel = open("noveltytests16.csv", "r")
frandom = open("randomtests16.csv", "r")
fmodnovel = open("modifiednoveltytests16.csv", "r")

results = []
results.append(readFile(faco))
results.append(readFile(fnovel))
results.append(readFile(frandom))
results.append(readFile(fmodnovel))


dfACO = pd.DataFrame.from_records([r.to_dict() for r in results[0]]).sort_values(by = ["tsp", "stag"])
dfNovel = pd.DataFrame.from_records([r.to_dict() for r in results[1]]).sort_values(by = ["tsp", "stag"])
dfRandom = pd.DataFrame.from_records([r.to_dict() for r in results[2]]).sort_values(by = ["tsp", "stag"])
dfModNovel = pd.DataFrame.from_records([r.to_dict() for r in results[3]]).sort_values(by = ["tsp", "stag"])

pcbDFACO = dfACO.loc[dfACO['tsp'] == "pcb442.tsp"]
pcbDFNovel = dfNovel.loc[dfNovel['tsp'] == "pcb442.tsp"]
pcbDFRandom = dfRandom.loc[dfRandom['tsp'] == "pcb442.tsp"]
pcbDFNovelMod = dfModNovel.loc[dfModNovel['tsp'] == "pcb442.tsp"]
#ratDF = dfACO.loc[dfACO['tsp'] == "rat783.tsp"]
#prDF = dfACO.loc[dfACO['tsp'] == "pr1002.tsp"]

pcbMeanACO = (pcbDFACO.groupby(['stag']).mean()).reset_index()
pcbMeanNovel = (pcbDFNovel.groupby(['stag']).mean()).reset_index()
pcbMeanRandom = (pcbDFRandom.groupby(['stag']).mean()).reset_index()
pcbMeanNovelMod = (pcbDFNovelMod.groupby(['stag']).mean()).reset_index()
#ratMean = (ratDF.groupby(['stag']).mean()).reset_index()
#prMean = (prDF.groupby(['stag']).mean()).reset_index()

#print(pcbMean.columns)
ax = pcbMeanACO.plot(x = 'stag', y = 'qual', marker = "o", title='pcb442 Results')
pcbMeanNovel.plot(x = 'stag', y = 'qual', marker = "o", ax = ax)
pcbMeanRandom.plot(x = 'stag', y = 'qual', marker = "o", ax = ax)
pcbMeanNovelMod.plot(x = 'stag', y = 'qual', marker = "o", ax = ax)

ax.set_xlabel("Stagnation Value")
ax.set_ylabel("Tour Length")
ax.legend(["ACO", "Novelty", "Random", "Modified Novel"]);
#ratMean.plot(x = 'stag', y = 'qual', marker = "o")
#prMean.plot(x = 'stag', y = 'qual', marker = "o")