import helper
import os
import time
import numpy as np
import cmath
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = np.argmax(energies > 10)
print(maxEnergyIndex)
fourierLength = int(param[fNum]["steps"] / 2)
sideLength = int(np.sqrt(param[fNum]["atoms"]))




def scatter(q, position):
    I_aa = np.zeros((3, fourierLength), dtype=complex)
    ic_sum = np.zeros(3, dtype=complex)
    for i in range(param[fNum]["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))

        I_aa[0, :] += np.fft.fft(q_dot_lattice * x[fNum, i, :].T)[:fourierLength]
        I_aa[1, :] += np.fft.fft(q_dot_lattice * y[fNum, i, :].T)[:fourierLength]
        I_aa[2, :] += np.fft.fft(q_dot_lattice * z[fNum, i, :].T)[:fourierLength]

        ic_sum[0] += cmath.exp(-1j * np.dot(q, position[i])) * x[fNum, i, 0]
        ic_sum[1] += cmath.exp(-1j * np.dot(q, position[i])) * y[fNum, i, 0]
        ic_sum[2] += cmath.exp(-1j * np.dot(q, position[i])) * z[fNum, i, 0]

    for i in range(3):
        I_aa[i, :] *= ic_sum[i]
        I_aa[i, :] = np.power(np.abs(I_aa[i, :]), 2)

    return I_aa.real

def runTransform():
    lattice_position = []
    for i in range(sideLength):
        for j in range(sideLength):
            lattice_position.append([j, i, 0])
    lattice_position = np.array(lattice_position).reshape(param[fNum]["atoms"], 3)

    I_total = np.zeros((sideLength, sideLength, maxEnergyIndex))

    q_x = np.array([0, 0, 0])
    eIndex = np.array([], dtype=np.float32)
    for i in range(sideLength):
        print(f"q_x: {i}")
        for j in range(sideLength):
            q_x = np.vstack((q_x, np.array([np.pi * (2 * i / sideLength - 1), np.pi * (2 * j / sideLength - 1), 0])))
            I_aa = scatter(q_x[-1, :], lattice_position)
            eIndex = np.append(eIndex, np.argmax(I_aa[0, :]))
            I_total[i, j, :] = I_aa[0, :maxEnergyIndex]

    try:
        os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
        os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat")
    except:
        print("Error while deleting files")
    I_total.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    eIndex.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat")

    return I_total, eIndex

recalculate = False
if recalculate:
    start = time.time()
    I_total, eIndex = runTransform()
    print(f'duration: {time.time()-start}')
else:
    I_total = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    I_total = I_total.reshape((sideLength, sideLength, maxEnergyIndex))
    eIndex = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat").astype(np.int32)

eIndex = eIndex[energies[eIndex].argsort()]
values, bins, _ = plt.hist(energies[eIndex], bins=144)
print(energies[eIndex])
eIndex = eIndex[np.where(values>=1)]

extent = [-np.pi, np.pi, -np.pi, np.pi]
fig, ax = plt.subplots(4, 4)
for i in range(4):
    for j in range(4):
        ax[i][j].imshow(I_total[:, :, int(eIndex[3 * i + j])], extent=extent)
        ax[i][j].set_aspect(abs(extent[1] - extent[0]) / abs(extent[3] - extent[2]))
        ax[i][j].set_title(f'E = {energies[int(eIndex[3 * i + j])]:.3f}')
# im = plt.imshow(I_total[:,:,eIndex], extent = extent)

# plt.plot(q[1:,0], helper.constants["J_to_meV"] * 4 * param[fNum]["J"]*3.5*(1-np.cos(q[1:,0])), 'w')


plt.show()
