import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec


n_instances = 4;
optimum = [5757191,6543610,6810665,8171677]

means_paco = [13.5,14.04,14.48,16.71]
means_partial = [5.45,5.82,5.81,7.20]

monaLisaSQ = (5855626+5854436+5855493+5854516+5855024+5855107+5855083+5855020+5855517+5854605)/10
vanGoghSQ = (6659450+6659519+6660497+6659129+6659261)/5
venusSQ = (6928149+6928519+6929370+6928488+6928252)/5
earringSQ = (8338885+8339629+8338089+8338282+8337626)/5

means_rpm =[monaLisaSQ,vanGoghSQ,venusSQ,earringSQ]
modified_means = []

for val, x in enumerate(means_rpm):
    diff = x-optimum[val]
    x = (diff/optimum[val])*100
    modified_means.append(x)
    


print(modified_means)


    
fig = plt.figure(figsize=(8,6))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
ax = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
plt.subplots_adjust(wspace = 0.02)
plt.xticks(fontsize = 8)

ax.grid(True,axis = "y", zorder = 0)
index = np.arange(n_instances)
bar_width = 0.3
opacity = 1

rects1 = ax.bar((index-bar_width/2), means_paco, bar_width,
                alpha=opacity, color='y', edgecolor='black', hatch = "///",zorder=3,
                label='P-ACO')

rects2 = ax.bar(index + bar_width/2, means_partial, bar_width,
                alpha=opacity, color='r',edgecolor='black', hatch = "|||",zorder=3,
                label='PartialACO')

rects3 = ax.bar(index + (1.5*bar_width), modified_means, bar_width,
                alpha=opacity, color='c',edgecolor='black',  hatch = "...",zorder=3,
                label='Restricted Pheromone Matrix')


ax.set_xlabel('Instance', fontsize=12)
ax.set_ylabel('% from shortest known')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(('monalisa100k', 'vangogh120k', 'venus140k', 'earring200k'))
ax.legend(loc=2, fontsize = 'small')

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

bplot1 = ax1.boxplot([modifiedData,modifiedData2], labels = ["Pheromone Map", "Heuristic"], widths=[0.7,0.7], zorder = 3, patch_artist = True, showmeans = True, meanline = True)


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
        
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
ax1.set_xlabel("Fallback Method")
ax1.set_ylabel("Percentage difference from shortest tour")
        
ax1.grid(True, axis = 'y', zorder= 0)
#plt.tight_layout()
plt.savefig("monalisacomparesmall.png",dpi=50)

plt.show()