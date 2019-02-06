import math
import matplotlib.pyplot as plt

def stats( data, num ):
	# (return min, max, mean, std deviation)
	tot = 0.0
	count = 0
	minVal = 1e30
	maxVal = -1e30

	n = len(data)
	if num != 0 and num < n:
		n = num

	for i in range(n):
		d = data[i]
		tot += d
		count += 1
		if d > maxVal:
			maxVal = d
		if d < minVal:
			minVal = d

	mean = tot / count
	var = 0.0
	for i in range(n):
		d = data[i]
		diff = d - mean
		diff *= diff
		var += diff
	var /= count
	sd = math.sqrt( var )

	return ( minVal, maxVal, mean, sd )

problems = ['lin318','pcb442','rat783','pr1002','pr2392']
optimum = [42029,50778,8806,259045,378032]

dataV = []
dataI = []

data = dataI
for algo in ['_ir_nn', '_vr_nn' ]:
	for p in problems:
		probData = []
		for suffix in ['20','40','70','100','0']:
			fileName = 'prodruns/'+p+algo+suffix+'.dat'
			print fileName
			irfile = open( fileName )
		
			solvals = []
			times = []
			for line in irfile:
				vals = line.split()
				sol = float(vals[0])
				ttour = float(vals[1])
				tpher = float(vals[2])
				solvals.append(sol)
				times.append(ttour+tpher)
			#calculate the stats
			sol_stats = stats( solvals, 0 )
			time_stats = stats( times, 0 )
			datapoint = [s for s in sol_stats]
			datapoint.append( time_stats[2] )
			probData.append( datapoint )
			
		data.append( probData )
	data = dataV

# put together some plots
numProblems = len(problems)

ord = [318.0,442.0,783.0,1002.0,2392.0] 

tiradotimes = [7.55,11.61,34.55,63.79,612.21]
reftimes = [0.0243, 0.05479,0.1952, 0.3307, 2.033]
reford = [318.0,442.0,783.0,1002.0, 2392.0]
tiradoplot = [t/256.0 for t in tiradotimes ]
colors = ['black', 'red']
plt.figure(num=None, figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')

#plt.text( 1, 1.75, 'vRoulette-2', color='black' )
#plt.text( 1, 1.70, 'vRoulette-1', color='red' )

for c in colors:
	if c == 'black':
		data = dataV
	else:
		data = dataI

	timePerIter = []
	for i in range( 0, numProblems ):
		timePerIter.append(data[i][0][4]/512.0)
#			err.append(v[3]/float(optimum[i]))
	plt.plot( ord, timePerIter, color=c, linewidth=2.0 )
	plt.scatter( ord, timePerIter, color=c )
#		deltay = -0.005
#		if ( c == 'red' and i == 4):
#			deltay -= 0.015
#		plt.text(  4.25, solq[4]+deltay, problems[i], color=c )
#		plt.xticks( ord, ['20','40','70','100','n'] )
plt.plot( ord, tiradoplot, color='green', linewidth=2.0)
plt.scatter( ord, tiradoplot, color='green')
plt.plot( reford, reftimes, color='blue', linewidth=2.0)
plt.scatter( reford, reftimes, color='blue')
ax=plt.gca()
ax.tick_params(axis='x',which='minor',bottom='off',top='off')
ax.set_yscale('log')
ax.set_xscale('log')

plt.xticks( ord, ['318','442','783','1002','2392'])

plt.text( 300, 0.015, 'CPU', color='blue')
plt.text( 1000, 0.7, 'Tirado', color='green')
plt.text( 300, 0.006, 'vRoulette-2', color='black')
plt.text( 1000, 0.02, 'vRoulette-1', color='red')
plt.ylim( 0.001, 4 )
plt.xlim( 200, 3000 )

plt.xlabel( 'Instance Size' )
plt.ylabel( 'Time per iteration (s)' )

plt.savefig('times.png',dpi=300)
plt.show()
