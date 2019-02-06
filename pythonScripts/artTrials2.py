import subprocess
import csv 
import random
import timeit

trials_per_seed = 1
num_seeds = 10
num_tsp = 1
seeds = []
tsp = ['mona-lisa100k.tsp']
for i in range(num_seeds):
    seeds.append(str(random.randint(0,99999)))
    
ant_commands = ['//home//staff//joshuap//Code//ACOCode//ACOCode//antmain', 'mona-lisa100k.tsp', '1000', '32', '40', '12345']

ant_commands_map = ['//home//staff//joshuap//Code//ACOCode//ACOCode//antmainmap', 'mona-lisa100k.tsp', '1000', '32', '40', '12345']

'''
with open('monalisastats2.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands[1] = "mona-lisa100k.tsp"
		randSeed = seeds[i]
		ant_commands[5] = randSeed
		output = subprocess.check_output( ant_commands)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])

with open('vangoughstats.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands[1] = "vangogh120K.tsp"
		randSeed = seeds[i]
		ant_commands[5] = randSeed
		output = subprocess.check_output( ant_commands)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])

with open('venusstats.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands[1] = "venus140K.tsp"
		randSeed = seeds[i]
		ant_commands[5] = randSeed
		output = subprocess.check_output( ant_commands)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])

with open('earringstats.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands[1] = "earring200K.tsp"
		randSeed = seeds[i]
		ant_commands[5] = randSeed
		output = subprocess.check_output( ant_commands)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])
'''

with open('monalisamapFINAL.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands_map[1] = "mona-lisa100k.tsp"
		randSeed = seeds[i]
		ant_commands_map[5] = randSeed
		output = subprocess.check_output( ant_commands_map)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(splitVals) == 6:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]) + " "+ str(splitVals[5]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]) + " " + str(splitVals[5])])
			elif len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])
'''
with open('vangoughstatsmap.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands_map[1] = "vangogh120K.tsp"
		randSeed = seeds[i]
		ant_commands_map[5] = randSeed
		output = subprocess.check_output( ant_commands_map)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])

with open('venusstatsmap.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands_map[1] = "venus140K.tsp"
		randSeed = seeds[i]
		ant_commands_map[5] = randSeed
		output = subprocess.check_output( ant_commands_map)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])

with open('earringstatsmap.csv', 'w') as csvfile:
	resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
	for i in range(num_seeds):
		resultWriter.writerow(['Run: ' + str(i)])
		ant_commands_map[1] = "earring200K.tsp"
		randSeed = seeds[i]
		ant_commands_map[5] = randSeed
		output = subprocess.check_output( ant_commands_map)
		vals = [s.strip() for s in output.splitlines()]
		for k in vals:
			splitVals = k.split()
			if len(k) > 0:
				print(str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4]))
				resultWriter.writerow([str(splitVals[0]) + " " + str(splitVals[1]) + " " + str(splitVals[2]) + " " +str(splitVals[3]) + " " + str(splitVals[4])])
'''
		
