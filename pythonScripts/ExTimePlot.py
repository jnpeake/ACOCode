import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(8,8))  # an empty figure with no axes

x = [318,442,783,1002,1577,2392,3795,5934,7397,11849]
plt.plot(x, [0.018191, 0.029170, 0.068259, 0.094123,0.176651, 0.371015, 0.785745, 2.088878, 3.388348,10.578780],'g',marker='D',ms='8',linewidth=2.5)
plt.plot(x, [0.0007342981054687499, 0.00136979412109375, 0.005648469277343749, 0.00901293005859375,0.020853308769531245, 0.046445145898437505, 0.14343100970703124, 0.4268420608593751, 0.7243076586914063,1.9750129865624997],'#339FFF',marker='s',ms='8',linewidth=2.5)
plt.plot(x, [0.0075768335, 1.4842791799999995/120, 2.9095891400000005/120, 3.6690928600000006/120,5.37324886/120, 12.2799055/120, 0.03991294222656249, 0.07261729623046877, 0.12465853121093749,0.27151259859374993],'y',marker='^',ms='8',linewidth=2.5)
plt.plot(x, [1.07299944/1024, 2.01855048/1024, 4.6500314199999995/1024, 3.0604545400000003/1024,4.8438548799999985/120, 6.851181539999999/120, 0.03991294222656249, 0.07261729623046877, 0.12465853121093749,0.27151259859374993],'b',marker='^',ms='8',linewidth=2.5)

plt.plot(x, [0.00052122369140625, 0.0007496182812500001, 0.0013743394726562499, 0.0018510785351562504 ,0.0035014900781250003, 0.007016489804687498, 0.013539661191406252, 0.027850285273437497, 0.04304569882812501,0.0974729433984375],'r',marker='o',ms='8',linewidth=2.5)

plt.legend(['ACOTSP - CPU','vRoulette-1 - Xeon Phi', 'VCSS - Xeon AVX2','VCSS - Xeon Phi AVX512'], loc='upper left',prop={'size': 16})
plt.xlabel('Domain Size')
plt.ylabel('Seconds per Iteration')
plt.yscale('log')
plt.xscale('symlog', linthreshy=0.001)
my_xticks = ['318', '442', '783', '1002','1577','2392','3795','5934','7397','11849']
plt.xticks(x, my_xticks)
ax = plt.gca()
plt.grid(True, axis='y',linestyle='--')

ax.tick_params(direction='in', length=6, width=2)
plt.savefig("ex_time.png",dpi=300)
plt.show()