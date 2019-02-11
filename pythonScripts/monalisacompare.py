#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:35:52 2019

@author: 13117809
"""

import matplotlib.pyplot as plt
import numpy as np

def get_stats( arr ):
    num = len(arr)
    mean = 0.0
    min = 1e20
    max = -1e20
    for v in arr:
        mean += v
        if v > max:
            max = v
        if v < min:
            min = v
            
    mean /= float(num)
    lq = sorted(arr)[num//4]
    median = sorted(arr)[num//2]
    uq = sorted(arr)[3*num//4]
    print ("Num: ",num, " Min: ",min, " LQ: " ,lq," Median: ",median, " UQ: ",uq," Max: ",max," Mean: ",mean)
    return (min, lq, median, uq, max, mean)



fig = plt.figure(figsize=(4,4))

optimum = 5757191

data = [5854279,5854685,5854125,5854225,5855313,5855764,5855340,5854885,5855490,5855301]
data2 = [5855626,5854436,5855493,5854516,5855024,5855107,5855083,5855020,5855517,5854605]

modifiedData = []
modifiedData2 = []


for x in data:
    diff = x-optimum
    x = (diff/optimum)*100
    modifiedData.append(x)
    
for x in data2:
    diff = x-optimum
    x = (diff/optimum)*100
    modifiedData2.append(x)


get_stats(modifiedData)
get_stats(modifiedData2)

bplot1 = plt.boxplot([modifiedData,modifiedData2], labels = ["Pheromone Map", "Heuristic"], widths=[0.7,0.7], zorder = 3, patch_artist = True, showmeans = True, meanline = True)


for val,patch in enumerate(bplot1['boxes']):
        patch.set_facecolor('white')
        patch.set_edgecolor('black')
        patch.set_linewidth(2)
        
for val,patch in enumerate(bplot1['medians']):
        patch.set_color('black')
        patch.set_linewidth(2)
        
for val,patch in enumerate(bplot1['means']):
        patch.set_color('red')
        patch.set_linewidth(2)
        
for val,patch in enumerate(bplot1['whiskers']):
        patch.set_linewidth(2)
        
for val,patch in enumerate(bplot1['caps']):
        patch.set_linewidth(2)
        
plt.xlabel("Fallback Method")
plt.ylabel("Percentage difference from shortest tour")
        
plt.grid(True, axis = 'y', zorder= 0)
plt.tight_layout()
plt.savefig("monalisacomparesmall.png",dpi=50)
