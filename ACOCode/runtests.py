import subprocess
import csv 
import random
import timeit


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
    perAnt = mean/10000
    print (num, min,lq,median,uq,max,mean, perAnt)
    return (min, lq, median, uq, max, mean, perAnt )

trials_per_seed = 1
num_seeds = 5
num_tsp = 1
seeds = []
#seeds = ['12345','54321','11111','22222','10101','99999','06994','98765','22222','33333']
#tsp = ['lin318.tsp','pcb442.tsp','rat783.tsp','pr1002.tsp','fl1577.tsp','pr2392.tsp','fl3795.tsp','rl5934.tsp','pla7397.tsp','rl11849.tsp']
#tsp = ['lin318.tsp','pcb442.tsp','rat783.tsp','pr1002.tsp','fl1577.tsp','pr2392.tsp']
#tsp = ['fl3795.tsp','rl5934.tsp','pla7397.tsp','rl11849.tsp']
#tsp = ['pla7397.tsp','rl11849.tsp']
#tsp = ['rl11849.tsp']
#tsp = ['rat783.tsp']
tsp = ['mona-lisa100k.tsp']
for i in range(num_seeds):
    seeds.append(str(random.randint(0,99999)))
    
ant_commands = ['//home//staff//joshuap//Code//ACOCode//ACOCode//antmain', 'a280.tsp', '10000', '32', '120', '12345']



totalTimeList = [];
totalLengthList = [];

with open('monalisastats.csv', 'w') as csvfile:
    resultWriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    resultWriter.writerow(['NNList', 'Seed'])
    resultWriter.writerow(['New'])
    print("NEW")
    for i in range(trials_per_seed):
        resultWriter.writerow([ant_commands[5]])
        for j in range (num_tsp):
            timeList = []
            totalTimeList = []
            lengthList = []
            ant_commands[1] = tsp[j]
            resultWriter.writerow(tsp[j])
            for k in range (num_seeds):
                randSeed = seeds[k]
                ant_commands[5] = randSeed
                start_time = timeit.default_timer()
                output = subprocess.check_output( ant_commands)
                vals = output.split()
                time = float(vals[-1])
                length = float(vals[-3])
                timeList.append(time)
                lengthList.append(length)
                end = timeit.default_timer() - start_time
                totalTimeList.append(end)
            resultWriter.writerow(timeList)
            resultWriter.writerow(lengthList)
            resultWriter.writerow(totalTimeList)
            resultWriter.writerow([str(get_stats(timeList))])
            resultWriter.writerow([str(get_stats(lengthList))])
            resultWriter.writerow([str(get_stats(totalTimeList))])
        totalTimeList.append((list(timeList), timeList[0]))
        totalLengthList.append((list(lengthList), lengthList[0]))
        resultWriter.writerow([ant_commands[5]])
        resultWriter.writerow(['--------------'])
        print("Seed:"+ant_commands[5])
        print("Times:"+(str(get_stats(timeList))))
        print("Lengths:"+(str(get_stats(lengthList))))
        
    ''' A comment
    resultWriter.writerow(['Original'])
    print("ORIGINAL")
    for i in range(trials_per_seed):
        resultWriter.writerow([ant_commands_original[5]])
        for j in range (num_tsp):
            timeList = []
            totalTimeList = []
            lengthList = []            
            ant_commands_original[1] = tsp[j]
            resultWriter.writerow(tsp[j])
            for k in range (num_seeds):
                randSeed = seeds[k]
                ant_commands_original[5] = randSeed
                start_time = timeit.default_timer()
                output = subprocess.check_output( ant_commands_original)
                vals = output.split()
                time = float(vals[-1])
                length = float(vals[-3])
                timeList.append(time)
                lengthList.append(length)
                end = timeit.default_timer() - start_time
                totalTimeList.append(end)
            resultWriter.writerow(timeList)
            resultWriter.writerow(lengthList)
            resultWriter.writerow(totalTimeList)
            resultWriter.writerow([str(get_stats(timeList))])
            resultWriter.writerow([str(get_stats(lengthList))])
            resultWriter.writerow([str(get_stats(totalTimeList))])
        totalTimeList.append((list(timeList), timeList[0]))
        totalLengthList.append((list(lengthList), lengthList[0]))
        resultWriter.writerow([ant_commands[5]])
        resultWriter.writerow(timeList)
        resultWriter.writerow(lengthList)
        resultWriter.writerow(['--------------'])
        print("Seed:"+ant_commands_original[5])
        print("Times:"+(str(get_stats(timeList))))
        print("Lengths:"+(str(get_stats(lengthList))))
'''
    


        


#print(vals)


