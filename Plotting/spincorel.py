import helper
import os, time
import numpy as np
import matplotlib.pyplot as plt
from turtle import shape
from matplotlib import colors
from scipy.signal import find_peaks
from matplotlib import rc

param, x, y, z = helper.getData()
fNum  = 0
nUnitCells = param[fNum]["nUnitCells"]

print(2*np.pi/param[fNum]["dt"])


def Intensity1(omega,starttimep = 0):
    Inten = 0j
    starttime = int(starttimep*param[fNum]["steps"])
    for n in range(param[fNum]["steps"]-1-starttime):
        correl01 = x[fNum, 0, n + starttime] * x[fNum, 1, n + 1 + starttime]
        Inten +=  np.exp(-1j * omega * n*param[fNum]["dt"]) * correl01 * param[fNum]["dt"]
    return Inten

def Intensity2(omega,starttimep = 0):
    Inten = 0j
    starttime = int(starttimep*param[fNum]["steps"])
    for n in range(param[fNum]["steps"]-1-starttime):
        correl01 = 0
        for i in range(nUnitCells):
            for j in range(nUnitCells):
                if j == nUnitCells - 1:
                    correl01 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, j*nUnitCells, n + 1 + starttime]
                else:
                    correl01 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, i + j*nUnitCells + 1, n + 1 + starttime]
        Inten +=  np.exp(-1j * omega * n*param[fNum]["dt"]) * correl01/nUnitCells**2 * param[fNum]["dt"]
    return Inten

def Intensity3(omega, starttimep=0):
    Inten = 0j
    starttime = int(starttimep*param[fNum]["steps"])
    for n in range(param[fNum]["steps"]-1-starttime):
        correl02 = x[fNum, 0, n + starttime] * x[fNum, nUnitCells, n + 1 + starttime]
        Inten +=  np.exp(-1j * omega * n*param[fNum]["dt"]) * correl02 * param[fNum]["dt"]
    return Inten

def Intensity4(omega,starttimep = 0):
    Inten = 0j
    starttime = int(starttimep*param[fNum]["steps"])
    for n in range(param[fNum]["steps"]-1-starttime):
        correl02 = 0
        for i in range(nUnitCells):
            for j in range(nUnitCells):
                if i == nUnitCells-1 and j == nUnitCells - 1:
                    correl02 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, 0, n + 1 + starttime]
                elif j == nUnitCells-1:
                    correl02 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, i + 1, n + 1 + starttime]
                elif n%2:
                    if i == nUnitCells-1:
                        correl02 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, (j+1)*nUnitCells, n + 1 + starttime]
                    else:
                        correl02 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, i + (j+1)*nUnitCells + 1, n + 1 + starttime]
                else:
                    correl02 += x[fNum, i + j*nUnitCells, n + starttime]*x[fNum, i + (j+1)*nUnitCells, n + 1 + starttime]
        Inten +=  np.exp(-1j * omega * n*param[fNum]["dt"]) * correl02/nUnitCells**2 * param[fNum]["dt"]
    return Inten
starts =0
stops = int(2*np.pi/param[fNum]["dt"]/5)

omegas = np.linspace(starts,stops,100)
#omegas = np.linspace(int(2*np.pi/param[fNum]["dt"]*5/11), int(2*np.pi/param[fNum]["dt"]*6/11),40)


int01 = np.array([])
int02 = np.array([])
start = time.time()
print('Calculating I:')
for i in omegas:
    if i < stops-1:
        print(f'progress: {((i-starts)/(stops-starts))*100:.2f} %', end='\r')
    else:
        print(f'progress: finshed')
    int01 = np.append(int01, Intensity2(i,0.5))
    int02 = np.append(int02, Intensity4(i,0.5))
print(f'duration: {time.time()-start}')

plt.figure()
plt.plot(omegas, np.absolute(int01),label="Abs[I01]")
plt.plot(omegas, np.absolute(int02),label="Abs[I02]")
plt.legend()
plt.show()