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

ord = range( 5 )

colors = ['black', 'red']
plt.figure(num=None, figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')

plt.text( 0, 1.28, 'vRoulette-2', color='black' )
plt.text( 0, 1.25, 'vRoulette-1', color='red' )

offsets = [1.02, 1.03, 1.04, 1.05, 1.07, ]
cpuqual = [1.034, 1.052, 1.080, 1.068, 1.17]
cpuord = [4.9,4.9,4.9,4.9,4.9]

plt.scatter( cpuord, cpuqual, color='blue' )
for i in range( 0, numProblems ):
	plt.text( 5.0, cpuqual[i], problems[i], color='blue')
plt.text( 4.8, 1.20, 'CPU', color='blue' )
for c in colors:
	if c == 'black':
		data = dataV
	else:
		data = dataI

	for i in range( 0, numProblems ):
		solq = []
		err = []
		for v in data[i]:
			solq.append(v[2]/float(optimum[i]))
			err.append(v[3]/float(optimum[i]))

		print len(ord), len(solq)
		plt.plot( ord, solq, color=c, linewidth=2.0 )
		plt.scatter( ord, solq, color=c )
		deltay = -0.005
		if ( c == 'red' and i == 4):
			deltay -= 0.015
		if ( c == 'red' ):
			plt.text(  -0.1, offsets[i], problems[i], color=c, horizontalalignment='right' )
		else:
			plt.text(  4.1,  solq[4]+deltay, problems[i], color=c )
		plt.xticks( range(6), ['20','40','70','100','n','CPU'] )

plt.xlim( -0.85, 5.9 )
plt.ylim( 1.0, 1.35 )
plt.xlabel( 'Nearest Neighbours' )
plt.ylabel( 'Solution Quality' )

plt.savefig('solqual.png',dpi=300)
plt.show()
