# Step number 0 calculation complete in 0.393974 sec with 999 vortices

import matplotlib.pyplot as plt
import re
from statistics import mean
from sys import argv

fName = ""
if len(argv) == 1:
	fName = "spawnRawOutput"
else:
	fName = argv[1]

f = open(fName, 'r')

totT = 0
step = 0

reg = re.compile("""^Step number (\d+) calculation complete in ([\d.]+)\D+(\d+)\D+$""")

# while True:
# 	line = f.readline()
# 	print(line)
# 	print(reg.search(line).group(1, 3))

matches = [reg.search(line) for line in f]
matches = [match for match in matches if match != None]
# print(matches)
data = {(float(x)*.01):float(y) for x, y in [match.group(1, 3) for match in matches]}

print("mean num vortices over the whole sim: " + str(mean(data.values())))

data = sorted(zip(data.keys(), data.values()), key = lambda datum: datum[0])

plt.plot([datum[0] for datum in data], [datum[1] for datum in data])


f.seek(0) #reset file

reg = re.compile("""spawning (\d+) vorts""")
matches = [reg.search(line) for line in f]
matches = [match for match in matches if match != None]
spawnRateData = [int(match.group(1)) for match in matches]

print("mean num spawns/step: " + str(mean(spawnRateData)))

f.seek(0) #reset file

reg = re.compile("""totMerges: (\d+)""")
matches = [reg.search(line) for line in f]
matches = [match for match in matches if match != None]
mergeRateData = [int(match.group(1)) for match in matches]

print("mean num merges/step: " + str(mean(mergeRateData)))


plt.show()