import subprocess

print ('Hello')
tsp = ['lin318.tsp','pcb442.tsp','rat783.tsp','pr1002.tsp','fl1577.tsp','pr2392.tsp','fl3795.tsp','rl5934.tsp','pla7397.tsp','rl11849.tsp']

ant_commands = ['//home//joshua//ACOCode//ACOCode//antmain', 'a280.tsp', '256', '32', '256', '12345']

output = subprocess.check_output( ant_commands)

print (output)
