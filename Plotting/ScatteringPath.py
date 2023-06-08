import helper
import numpy as np
from scipy.signal import find_peaks
import os
import time
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = int(param[fNum]["steps"]//2+1)
fourierLength = int(param[fNum]["steps"] / 2)

pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
while os.path.exists(pathFourier):
    os.remove(pathFourier)
    fNum += 1
    pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
fNum = 0

def getTransform():
    sideLength = param[fNum]["nUnitCells"]
    latticePosition = helper.getPositions(fNum)
    scatterLength = int(sideLength/2)
    if param[fNum]["geometry"] == "line":
        qScatter = helper.scatterLine(scatterLength)

    elif param[fNum]["geometry"] == "square":
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/4)
        qScatter = helper.scatterSquare(scatterLength)

    elif param[fNum]["geometry"] == "triangle":
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/2)
        qScatter = helper.scatterTriangle(scatterLength)
    
    elif param[fNum]["geometry"] == "kagome":
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/2)
        qScatter = helper.scatterKagome(scatterLength)

    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
    if not os.path.exists(pathFourier):
        print("Calculating")
        start = time.time()
        print("x: ")
        I_total_x = helper.runTransform(latticePosition, x[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        print("y: ")
        I_total_y = helper.runTransform(latticePosition, y[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        print("z: ")
        # I_total_z = helper.runTransform(latticePosition, z[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        I_total = I_total_x**2 + I_total_y**2# + I_total_z**2
        I_total.tofile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
        print(f'duration: {time.time()-start}')
    else:
        print("Pulling fourier data from file")
        I_total = np.fromfile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
    
    return I_total, nScatter, sideLength

fig, ax = plt.subplots(2)
point = 1
usePeaks = True
while fNum < x.shape[0]:
    I_total, nScatter, sideLength = getTransform()
    # xData = np.linspace(0,nScatter-1,nScatter)
    xData, yData, zData = [], [], []
    maxF = np.log10(np.max(I_total))
    for i in range(nScatter):
        if usePeaks:
            peaks, _ = find_peaks(I_total[i, :], height=np.max(I_total[i, :])/1000, distance=50)
            # Find peaks cannot check if the peaks at index 0
            if I_total[i, 0] >= np.max(I_total[i, :])/10:
                yData.append(np.abs(energies[0]))
                xData.append(i)
                zData.append(max([np.log10(np.max(I_total[i,0]))-maxF, -5]))

            # If there are no peaks take the highest value
            if peaks.size == 0:
                yData.append(np.abs(energies[np.argmax(I_total[i,:])]))
                xData.append(i)
                zData.append(max([np.log10(np.max(I_total[i,:]))-maxF, -5]))

            # Plot all the peaks
            else:
                for peak in peaks:
                    yData.append(np.abs(energies[peak]))
                    xData.append(i)
                    zData.append(max([np.log10(np.max(I_total[i, peak]))-maxF, -5]))
        
        else:
            yData.append(np.abs(energies[np.argmax(I_total[i,:])]))
            xData.append(i)
            zData.append(max([np.log10(np.max(I_total[i,:]))-maxF, -5]))

    # ax[0].plot(xData,yData, helper.constants["colors"][fNum]+'o')
    helper.plotTheory(param[fNum]["geometry"], sideLength, param[fNum], ax[0], fNum)
    im = ax[0].scatter(xData,yData,c=zData,cmap="Wistia", zorder=10)
    fig.colorbar(im, ax=ax[0])

    helper.plotPeaks(ax[1], energies[:maxEnergyIndex-1], I_total[point,:])
    ax[1].plot(energies[:maxEnergyIndex-1], I_total[point,:maxEnergyIndex-1], c= helper.constants["colors"][fNum])
    helper.plotSmooth(ax[1], energies[:maxEnergyIndex-1], I_total[point,:])
    fNum = fNum + 1

fNum = fNum - 1

ax[1].set_yscale("log")
ax[0].axvline(point, c='r')
ax[0].set_ylabel('Energy [meV]')
ax[0].legend(loc="lower right")

if param[fNum]["J"] > 0:
    ax[0].set_title('Ferromagnetic spin chain')
else:
    ax[0].set_title('Antiferromagnetic spin chain')
helper.plotMarker(param[fNum]["geometry"], sideLength, ax[0])

plt.show()
