import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(8,8))  # an empty figure with no axes

x = [318,442,783,1002,1577,2392,3795,5934,7397,11849]
plt.plot(x, [0.018191, 0.029170, 0.068259, 0.094123,0.176651, 0.371015, 0.785745, 2.088878, 3.388348,10.578780],'g',marker='o',linewidth=2.5)
plt.plot(x, [0.0002829267187500001, 0.000433415234375, 0.0015416902148437505, 0.00245905220703125,0.005466861738281249, 0.011835663554687499, 0.0356230121875, 0.1109166513671875, 0.1840920808984375,0.5032022111132812],'b',marker='o',linewidth=2.5)
plt.plot(x, [0.0006649209093302488, 0.00082676400390625, 0.0017667200390625, 0.0018471205859374997,0.003971166386718749, 0.007027976621093748, 0.017914368476562505, 0.044908897910156266, 0.04818707912109376,0.19480160562500004],'r',marker='o',linewidth=2.5)

plt.legend(['ACOTSP - CPU','vRoulette-1 - Xeon Phi', 'VCSS - Xeon Phi'], loc='upper left',prop={'size': 16})
plt.xlabel('Domain Size')
plt.ylabel('Seconds per Iteration')
plt.yscale('log')
plt.xscale('symlog', linthreshy=0.001)
my_xticks = ['318', '442', '783', '1002','1577','2392','3795','5934','7397','11849']
plt.xticks(x, my_xticks)
ax = plt.gca()

ax.tick_params(direction='in', length=6, width=2)
plt.savefig("execution_time.png",dpi=300)
plt.show()