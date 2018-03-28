import matplotlib.pyplot as plt
import csv

def plot_box( ax, x, boxcol, boxalpha, width, yvals ):
    box_x = [x-width,x-width,x+width,x+width,x-width]
    box_y = [yvals[1], yvals[3], yvals[3], yvals[1], yvals[1]]
    ax.fill( box_x, box_y, boxcol )
    ax.plot( box_x, box_y, c=boxcol )
    ax.plot( [x-width,x+width], [yvals[2],yvals[2]], c='black')
    ax.plot( [x,x], [yvals[1],yvals[0]], c=boxcol)
    ax.plot( [x,x], [yvals[3],yvals[4]], c=boxcol)
    ax.plot( [x-width,x+width],[yvals[0],yvals[0]],c=boxcol)
    ax.plot( [x-width,x+width],[yvals[4],yvals[4]],c=boxcol)

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
    return (min, lq, median, uq, max, mean )

# optimum tour lengths 
optimum = { 'lin318' : 42029.0, 'pcb442' : 50778.0, 'rat783' : 8806.0, 'pr1002' : 259045.0, 'fl1577' : 22249.0,
			'pr2392' : 379032.0, 'fl3795' : 28772.0, 'rl5934' : 556045.0, 'pla7397' : 23260728.0, 'rl11849' : 923288.0 }

instance_names = [v for v in optimum]

# CPU run stats
cpu_stats = {}

# solution quality - dictionary with instance name as key, list of tour lengths as value
sol_qual_no_nn = {}
sol_qual_nn = {}

# xeon phis stats
stats_no_nn = {}
stats_nn = {}

sol_qual = sol_qual_nn # results with nn list are first
cur_instance = None
countdown = 0


with open('test.csv') as csvfile:
	data_reader = csv.reader(csvfile)
	for line in data_reader:
		# have we hit the 'Original' data yet?
		if 'Original' in line:
			sol_qual = sol_qual_no_nn

		# try to identify this as a tsp file
		munged = ''.join(line)
		if '.tsp' in munged:
			# ok, we have an instance. We want the data in the line after next 
			cur_instance = munged.replace('.tsp','')
			countdown = 2

		if countdown == 0 and cur_instance != None:
			sol_qual[cur_instance] = [float(val) for val in line]
			cur_instance = None
		
		if countdown > 0: 
			countdown -= 1

# read CPU result stats
with open('poincare_runs.txt') as in_file:
	for line in in_file:
		vals = line.split()
		if vals[0] in instance_names: # look for the instance name in the first column
			cpu_stats[vals[0]] = [float(vals[i]) for i in [3,4,5,6,7,2]] # ignore first two, and last column


# calculate stats for xeon phi runs
for name in instance_names:
	stats_no_nn[name] = [v for v in get_stats(sol_qual_no_nn[name])]
	stats_nn[name] = [v for v in get_stats(sol_qual_nn[name])]

# normalise all stats to optimum
for stats in [cpu_stats, stats_no_nn, stats_nn]:
	for name in instance_names:
		for i in range(len(stats[name])):
			stats[name][i] /= optimum[name]


#
# we have all the data, now plot...
#

plt.figure(num=None, figsize=(8,6), facecolor='w', edgecolor='k' )
ax = plt.gca()
ax.set_xticks([])
plt.xlim(0,10)
plt.ylim(1.0,1.5)


for i in range(0,10,2):
    x = float(i)
    ax.fill( [x,x,x+1,x+1,x],[1,1.5,1.5,1,1], color='gray', alpha=0.2 )

x = 0
y = 1.3
for name in instance_names:
	plt.text( x+0.5, y, name, {'ha':'center', 'va':'bottom'}, rotation=90, fontsize=12)
	plot_box( ax, x+0.2, 'green', 1.0, 0.1, cpu_stats[name] )
	plot_box( ax, x+0.5, 'red', 1.0, 0.1, stats_no_nn[name] )
	plot_box( ax, x+0.8, 'blue', 1.0, 0.1, stats_nn[name] )
	x += 1
	if name == 'pr1002':
		y += 0.08 

# legend
lx = 0.8
ly = 1.38
lw = 3
lh = 0.1

ax.fill( [lx, lx, lx+lw, lx+lw, lx],[ly,ly+lh,ly+lh,ly,ly], color='white' )
ax.plot( [lx, lx, lx+lw, lx+lw, lx],[ly,ly+lh,ly+lh,ly,ly], color='black' )
plt.text( lx + lw*0.1, ly + lh * 0.8, 'CPU', color='green')
plt.text( lx + lw*0.1, ly + lh * 0.5, 'Xeon Phi, no NN', color='red')
plt.text( lx + lw*0.1, ly + lh * 0.2, 'Xeon Phi, with NN', color='blue')
plt.xlabel('Instance',fontsize=14)
plt.ylabel('Solution Quality',fontsize=14)
plt.savefig("solution_quality.png",dpi=300)
plt.show()
