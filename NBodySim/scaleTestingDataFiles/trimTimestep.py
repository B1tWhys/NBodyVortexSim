# usage:
# python trimTimestep.py newDomainSize originFile outputFile

from sys import argv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.image#  as image
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
    def __init__(self, xpos, ypos, xVel, yVel, tVel):
        self.xpos = xpos
        self.ypos = ypos
        self.xVel = xVel
        self.yVel = yVel
        self.tVel = tVel # total vel. magnitude

class Vortex(PointObj):
    def __init__(self, index, xpos, ypos, xVel, yVel, tVel, vorticity, spawnStep):
        super(Vortex, self).__init__(xpos, ypos, xVel, yVel, tVel)
        self.vorticity = vorticity
        self.index = index
        self.spawnStep = spawnStep

class Tracer(PointObj):
    def __init__(self, index, xpos, ypos, xVel, yVel, tVel):
        super(Tracer, self).__init__(xpos, ypos, xVel, yVel, tVel)
        self.index = index

inF = open(argv[2], 'r')
fSize = os.path.getsize(argv[2])

buffer = ""

char = inF.read(1)

cutoff = float(argv[1])

assert (char == '\x1D'), "Missing group seperator (\x1D) at byte %i"%(inF.tell())

line = inF.readline().strip()
strArr = line.split(',')
ts = Timestep(int(strArr[0]), float(strArr[1]), int(strArr[2]))
numV = int(strArr[3])
numT = int(strArr[4])
# timesteps.append(ts)
sep = inF.read(1)

assert (sep == '\x1E'), "Missing record seperator (\x1E) at byte %i. Value is %s instead"%(inF.tell(), sep)
for i in range(numV):
    line = inF.readline().strip()
    strArr = line.split(',')
    index = strArr.pop(0)
    xpos = float(strArr.pop(0))
    ypos = float(strArr.pop(0))
    xVel = float(strArr.pop(0))
    yVel = float(strArr.pop(0))
    tVel = 0.0
    if (len(strArr) == 8):
        tVel = float(strArr.pop(0))
    vorticity = float(strArr.pop(0))
    spawnStep = int(strArr.pop(0))
    
    vort = Vortex(index, xpos, ypos, xVel, yVel, tVel, vorticity, spawnStep)
    ts.vorts.append(vort)

assert (inF.read(1) == '\x1E'), "Missing record seperator (\x1E) at byte %i"%(inF.tell())
for i in range(numT):
    line = inF.readline().strip()
    strArr = line.split(',')
    
    index = int(strArr.pop(0))
    xpos = float(strArr.pop(0))
    ypos = float(strArr.pop(0))
    xVel = float(strArr.pop(0))
    yVel = float(strArr.pop(0))
    tVel = 0.0
    if (len(strArr) == 8):
        tVel = float(strArr.pop(0))    
    
    tracer = Tracer(index, xpos, ypos, xVel, yVel, tVel)
    ts.tracers.append(tracer)

print "cutting off at %f"%(cutoff)
for i in range(len(ts.vorts)-1, -1, -1):
    vort = ts.vorts[i]
    if (vort.xpos > cutoff):
        ts.vorts.remove(vort)
    elif (vort.ypos > cutoff):
        ts.vorts.remove(vort)

for i, vort in enumerate(ts.vorts):
    vort.index = i

outF = open(argv[3], 'w')
outF.write("\x1D0,0.0,%i,%i,%i\n\x1E"%(ts.seedVal, len(ts.vorts), len(ts.tracers)))
for vort in ts.vorts:
    outF.write("%i,%f,%f,%f,%f,%f,%i\n"%(vort.index, vort.xpos, vort.ypos, vort.xVel, vort.yVel, vort.vorticity, vort.spawnStep))
outF.write("\x1E")
for tracer in ts.tracers:
    outF.write("%i,%f,%f,%f,%f\n"%(tracer.index, tracer.xpos, tracer.ypos, tracer.xVel, tracer.yVel))
outF.write("\x1D")
inF.close()
outF.close()
