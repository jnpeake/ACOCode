import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(8,8))  # an empty figure with no axes

x = [318,442,783,1002,1577,2392]
plt.plot(x, [42863,53617,9599,308699,31201,581351],'g',marker='D',ms='8',linewidth=2.5)
plt.plot(x, [42763,53222,9398,287924,24307,453656],'#339FFF',marker='s',ms='8',linewidth=2.5)

plt.legend(['Default Value','Standard VCSS'], loc='upper left',prop={'size': 16})
plt.xlabel('Domain Size')
plt.ylabel('Seconds per Iteration')
plt.xscale('symlog', linthreshy=0.001)
my_xticks = ['318', '442', '783', '1002','1577','2392','3795','5934','7397','11849']
plt.xticks(x, my_xticks)
ax = plt.gca()
plt.grid(True, axis='y',linestyle='--')

ax.tick_params(direction='in', length=6, width=2)
plt.savefig("ex_time.png",dpi=300)
plt.show()