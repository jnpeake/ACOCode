import matplotlib.pyplot as plt
import numpy as np
import re

filenameML = "MonaLisaStats.csv"
filenameVG = "VanGoghStats.csv"
filenameV = "VenusStats.csv"
filenameE = "EarringStats.csv"
fileML = open(filenameML, "r")
fileVG = open(filenameVG, "r")
fileV = open(filenameV, "r")
fileE = open(filenameE, "r")
files = [fileML,fileVG,fileV,fileE]

optimum = [5757191,6543610,6810665,8171677]
averageLists = []
for i, lines in enumerate(files):
    
    averageListPre = np.zeros(101, dtype=float)
    averageList = np.zeros(101, dtype=float)
    for line in lines:
        line = re.sub('b', '', line)
        splitLine = line.split()
        if splitLine[1] != "'tour:'" and splitLine[0] != "Run:":
           iteration = splitLine[0]
           solqual = splitLine[1]
           iteration = re.sub('\'', '', iteration)
           solqual = re.sub('\'', '', solqual)
           iterationNum = int(int(iteration)/10)
           averageListPre[int(iterationNum)] +=  float(solqual)
           
        elif splitLine[1] == "'tour:'":
           iteration = '1000'
           solqual = splitLine[2]
           iteration = re.sub('\'', '', iteration)
           solqual = re.sub('\'', '', solqual)
           
           averageListPre[100] +=  float(solqual)
           
        for val, a in enumerate(averageListPre):
               if i == 0:
                   a = a/10
               else:
                   a = a/5
               diff = a-optimum[i]
               a = (diff/optimum[i])*100
               averageList[val] = a
           
    averageLists.append(averageList)
       

fig = plt.figure(figsize=(8,8))  # an empty figure with no axes

x = [0,100,200,300,400,500,600,700,800,900,1000]
ax = fig.add_subplot(111)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)



ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')



ax1.set_title('mona-lisa100k')
ax2.set_title('vangogh120k')
ax3.set_title('venus140k')
ax4.set_title('earring200k')
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
ax4.grid(True)
ax1.plot(x,[averageLists[0][0],averageLists[0][10],averageLists[0][20],averageLists[0][30],averageLists[0][40],averageLists[0][50],averageLists[0][60],averageLists[0][70],averageLists[0][80],averageLists[0][90],averageLists[0][100]], color='red')
ax2.plot(x,[averageLists[1][0],averageLists[1][10],averageLists[1][20],averageLists[1][30],averageLists[1][40],averageLists[1][50],averageLists[1][60],averageLists[1][70],averageLists[1][80],averageLists[1][90],averageLists[1][100]], color='red')
ax3.plot(x,[averageLists[2][0],averageLists[2][10],averageLists[2][20],averageLists[2][30],averageLists[2][40],averageLists[2][50],averageLists[2][60],averageLists[2][70],averageLists[2][80],averageLists[2][90],averageLists[2][100]], color='red')
ax4.plot(x,[averageLists[3][0],averageLists[3][10],averageLists[3][20],averageLists[3][30],averageLists[3][40],averageLists[3][50],averageLists[3][60],averageLists[3][70],averageLists[3][80],averageLists[3][90],averageLists[3][100]], color='red')

ax.set_ylabel("Percentage difference from shortest tour")
ax.set_xlabel("Iterations")

plt.savefig("convergencesmall.png",dpi=50)       