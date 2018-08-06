# usage:
# python dataAnalyzer.py path/to/rawDataFile [out/file/path]
# output file path defaults to ./data/plots_1.png

from sys import argv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.image as image
import matplotlib.pyplot as plt
# import array
from collections import deque
import statistics as stats
import math

class Timestep(object):
    def __init__(self, stepNum, time, seedVal):
        self.index = stepNum
        self.time = time
        self.seedVal = seedVal
        self.vorts = []
        self.tracers = []

class PointObj(object):
    def __init__(self, xPos, yPos, xVel, yVel, tVel):
        self.xpos = xPos
        self.ypos = yPos
        self.xVel = xVel
        self.yVel = yVel
        self.tVel = tVel # total vel. magnitude

class Vortex(PointObj):
    def __init__(self, index, xPos, yPos, xVel, yVel, tVel, vorticity, spawnStep):
        super(Vortex, self).__init__(xPos, yPos, xVel, yVel, tVel)
        self.vorticity = vorticity
        self.index = index
        self.spawnStep = spawnStep

class Tracer(PointObj):
    def __init__(self, index, xPos, yPos, xVel, yVel, tVel):
        super(Tracer, self).__init__(xPos, yPos, xVel, yVel, tVel)
        self.index = index

f = open(argv[1], 'r')
fSize = os.path.getsize(argv[1])

buffer = ""

# timesteps = []
times = deque()
numVorts = deque()
numPosVorts = deque()
numNegVorts = deque()
gamma_pos = deque()
gamma_neg = deque()
gamma_tot = deque()

while True:
    char = f.read(1)
    
    if (char == ""): break
    
    assert (char == '\x1D'), "Missing group seperator (\x1D) at byte %i"%(f.tell())
    
    line = f.readline().strip()
    strArr = line.split(',')
    ts = Timestep(int(strArr[0]), float(strArr[1]), int(strArr[2]))
    if (ts.index % 100 == 0):
        print("computing timestep: %i"%ts.index)
#    if ts.index == 10000:
#        break;
    numV = int(strArr[3])
    numT = int(strArr[4])
    # timesteps.append(ts)
    sep = f.read(1)
    
    assert (sep == '\x1E'), "Missing record seperator (\x1E) at byte %i. Value is %s instead"%(f.tell(), sep)
    for i in range(numV):
        line = f.readline().strip()
        strArr = line.split(',')
         
        index = int(strArr.pop(0))
        xPos = float(strArr.pop(0))
        yPos = float(strArr.pop(0))
        xVel = float(strArr.pop(0))
        yVel = float(strArr.pop(0))
        tVel = 0.0
        if (len(strArr) == 8):
            tVel = float(strArr.pop(0))
        vorticity = float(strArr.pop(0))
        spawnStep = int(strArr.pop(0))
        
        vort = Vortex(index, xPos, yPos, xVel, yVel, tVel, vorticity, spawnStep)
        ts.vorts.append(vort)
    
    assert (f.read(1) == '\x1E'), "Missing record seperator (\x1E) at byte %i"%(f.tell())
    for i in range(numT):
        line = f.readline().strip()
        strArr = line.split(',')
        
        index = int(strArr.pop(0))
        xPos = float(strArr.pop(0))
        yPos = float(strArr.pop(0))
        xVel = float(strArr.pop(0))
        yVel = float(strArr.pop(0))
        tVel = 0.0
        if (len(strArr) == 8):
            tVel = float(strArr.pop(0))    
        
        tracer = Tracer(index, xPos, yPos, xVel, yVel, tVel)
        ts.tracers.append(tracer)
    
    ts.vorts = sorted(ts.vorts, key = lambda vort: vort.index)
    ts.tracers = sorted(ts.tracers, key = lambda tracer: tracer.index)
    
    totVortCount = len(ts.vorts)
    posVortCount = len([vort for vort in ts.vorts if vort.vorticity > 0])
    negVortCount = totVortCount - posVortCount
    posGammaSum = sum([vort.vorticity for vort in ts.vorts if vort.vorticity > 0])/(math.pi*2)
    negGammaSum = sum([vort.vorticity for vort in ts.vorts if vort.vorticity < 0])/(math.pi*2)
    gammaSum = (posGammaSum + negGammaSum)/(math.pi*2)
    
    times.append(ts.index*.01)
    numVorts.append(totVortCount)
    numPosVorts.append(posVortCount)
    numNegVorts.append(negVortCount)
    gamma_pos.append(posGammaSum)
    gamma_neg.append(negGammaSum)
    gamma_tot.append(gammaSum)
#     print("calculated step %i"%ts.
#     print(len([vort.vorticity for vort in ts.vorts if vort.vorticity > .001 and vort.vorticity < .002]))
    
########### draw the graphs
avgN = stats.mean(numVorts)
print("average vort count: %f"%avgN)

# plt.yticks(list(range(0, 800, 50)))
fig, plots = plt.subplots(nrows=2)
ax1 = plots[0]
ax1.plot(times, numVorts, 'k-')
ax1.plot(times, numPosVorts, 'r-')
ax1.plot(times, numNegVorts, 'b-')
ax1.set(xlabel='time', ylabel='num vortices')
ax1.grid()
# plt.show()

# fig, ax = plt.subplots()
ax2 = plots[1]
# plt.yticks(np.arange(-800, 800, 50))
ax2.plot(times, gamma_tot, 'k-')
ax2.plot(times, gamma_pos, 'r-')
ax2.plot(times, gamma_neg, 'b-')
ax2.grid()

ax1.set(xlabel='time', ylabel='num vortices')
ax2.set(xlabel='time', ylabel='Gamma')

plotsfName = argv[2] if (len(argv) == 3) else "./data/plots_1.png"

plt.savefig(plotsfName, dpi=300)
# plt.show()
