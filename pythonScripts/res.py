import matplotlib.pyplot as plt

datfile = open("plot.dat")
names = [ 'd198', 'a280', 'lin318', 'pcb442', 'rat783', 'pr1002', 
		  'pcb1173', 'd1921', 'pr2392' ]
qir = []
sigir = []
qr = []
sigr = []

lim0 = 1.02
lim1 = 1.42

for line in datfile:
	tokens = line.split()
	fvals = [ float(x) for x in tokens ]
	qr.append( fvals[0] )
	sigr.append( fvals[2] )
	qir.append( fvals[3] )
	sigir.append( fvals[5] )

plt.plot( [lim0, lim1], [lim0, lim1], '--' )
ax = plt.gca()
ax.set_xticks( [1,1.1,1.2,1.3,1.4])
ax.set_yticks( [1,1.1,1.2,1.3,1.4])
ax.set_xlabel( "$Q_R$", fontsize=18 )
ax.set_ylabel( "$Q_{IR}$", fontsize=18 )

for i in range (0, len(names) ):
	text = names[i]
	x = qr[i]  
	y = qir[i] + 0.075
	plt.text( x, y, text, rotation=90, fontsize=15, bbox={'facecolor':'white', 'edgecolor':'white'} )

plt.errorbar( qr, qir, xerr=sigr, yerr=sigir, fmt='o' )
plt.show()


