import helper
import os
import time
import numpy as np
import cmath
import matplotlib.pyplot as plt

fNum = 0
param, x, y, z = helper.getData2(fNum)


frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies/np.pi**2
maxEnergyIndex = np.argmax(energies > 20)
#maxEnergyIndex = int(param[fNum]["steps"]//2+1)
print(maxEnergyIndex)
fourierLength = int(param[fNum]["steps"] / 2)
sideLength = int(np.sqrt(param[fNum]["atoms"]))

scalling1 = np.pi/3*8
scalling2 = np.pi/3*8*np.sqrt(3)/2
qmult = 1
resol = sideLength*2
maxq1 = scalling1*qmult
minq1 = -scalling1*qmult
maxq2 = scalling2*qmult
minq2 = -scalling2*qmult
recalculate = True


xtlab = np.linspace(-1*qmult, 1*qmult, 2*qmult+1,dtype=np.int32)
xtic = np.linspace(-scalling1, scalling1, np.size(xtlab))
ytlab = np.linspace(-2*qmult, 2*qmult, 4*qmult+1,dtype=np.int32)
ytic = np.linspace(-scalling2*qmult, scalling2*qmult, np.size(ytlab))

pathFourier = helper.getPath(helper.constants["pathFourier"], 0)
pathEnergyIndex = helper.getPath(helper.constants["pathEnergyIndex"], 0)


def getFasterTransform(sx,sy,sz,qScatter, i = 0):
    latticePosition = helper.getPositions(i)
    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], i)
    makenew = True #if i want to calculate new fourier
    if not os.path.exists(pathFourier) or makenew:
        print("Calculating")
        start = time.time()
        print("x: ")
        I_total_x = helper.getF_to_I_dp(sx[:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        print("y: ")
        I_total_y = helper.getF_to_I_dp(sy[:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        print("z: ")
        I_total_z = helper.getF_to_I_dp(sz[:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        I_total =   I_total_x**2 + I_total_y**2 + I_total_z**2
        I_total.tofile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
        print(f'total duration: {time.time()-start}')
    else:
        print("Pulling fourier data from file")
        I_total = np.fromfile(pathFourier)
        I_total = I_total.reshape((nScatter, -1))
    return I_total
'''
def scatter3(q, position):
    I_aa = np.zeros((3, fourierLength), dtype=complex)
    #ic_sum = np.zeros(3, dtype=complex)
    for i in range(param[fNum]["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))

        I_aa[0, :] += np.fft.fft(q_dot_lattice * x[fNum, i, :].T)[:fourierLength]
        I_aa[1, :] += np.fft.fft(q_dot_lattice * y[fNum, i, :].T)[:fourierLength]
        I_aa[2, :] += np.fft.fft(q_dot_lattice * z[fNum, i, :].T)[:fourierLength]

        #ic_sum[0] += cmath.exp(-1j * np.dot(q, position[i])) * x[fNum, i, 0]
        #ic_sum[1] += cmath.exp(-1j * np.dot(q, position[i])) * y[fNum, i, 0]
        #ic_sum[2] += cmath.exp(-1j * np.dot(q, position[i])) * z[fNum, i, 0]

    for i in range(3):
        #I_aa[i, :] *= ic_sum[i]
        I_aa[i, :] = np.power(abs(I_aa[i, :]), 2)

    return np.sqrt(I_aa.real[0]**2 + I_aa.real[1]**2 + I_aa.real[2]**2)
    
def scatter2(q, position):
    I_xx = helper.Fourier(x[fNum,:,:],param[fNum])
    I_yy = helper.Fourier(y[fNum,:,:],param[fNum])
    I_zz = helper.Fourier(z[fNum,:,:],param[fNum])
    
    for i in range(param[fNum]["atoms"]):
        I_xx[i,:] *= np.exp(1j * np.dot(q, position[i])) * x[fNum, i, 1000]
        I_yy[i,:] *= np.exp(1j * np.dot(q, position[i])) * y[fNum, i, 1000]
        I_zz[i,:] *= np.exp(1j * np.dot(q, position[i])) * z[fNum, i, 1000]
    """
    I_xx = np.power(I_xx, 2)
    I_xx = np.power(I_xx, 2)
    I_xx = np.power(I_xx, 2)
    I_xx = np.sum(I_xx, axis = 0)
    I_yy = np.sum(I_yy, axis = 0)
    I_zz = np.sum(I_zz, axis = 0)
    """
    I_xx = abs(I_xx[:maxEnergyIndex])
    I_yy = abs(I_yy[:maxEnergyIndex])
    I_zz = abs(I_zz[:maxEnergyIndex])
    return np.sqrt(I_xx**2 + I_yy**2 + I_zz**2)

def scatter(q, position):
    I_aa = np.zeros((3, fourierLength), dtype=complex)
    for i in range(param[fNum]["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))

        I_aa[0, :] += np.fft.fft(q_dot_lattice * x[fNum, i, :].T)[:fourierLength]
        I_aa[1, :] += np.fft.fft(q_dot_lattice * y[fNum, i, :].T)[:fourierLength]
        I_aa[2, :] += np.fft.fft(q_dot_lattice * z[fNum, i, :].T)[:fourierLength]

    q_norm = q/np.linalg.norm(q)
    return (1-q_norm[0]**2)*np.power(abs(I_aa[0]), 2)+(1-q_norm[1]**2)*np.power(abs(I_aa[1]), 2)+(1-q_norm[2]**2)*np.power(abs(I_aa[2]), 2)

def runTransform(minrec = -np.pi, maxrec = np.pi, resolution = sideLength):
    lattice_position = helper.getPositions(fNum)
    lattice_position = np.array(lattice_position).reshape(param[fNum]["atoms"], 3)

    I_total = np.zeros((resolution, resolution, maxEnergyIndex))

    q_x = np.array([0, 0, 0])
    eIndex = np.array([], dtype=np.int32)
    for i in range(resolution):
        for j in range(resolution):
            q_x = np.vstack((q_x, np.array([minrec+i/resolution*(maxrec-minrec),minrec+j/resolution*(maxrec-minrec),0])))
            I_aa = scatter(q_x[-1, :], lattice_position)
            eIndex = np.append(eIndex, np.argmax(I_aa[:]))
            I_total[i, j, :] = I_aa[:maxEnergyIndex]
            print(f'q_x {i} progress: {j/resolution*100:.2f} %. total progress {i/resolution*100:.2f} %', end='\r')
    
    I_total.tofile(pathFourier)
    eIndex.tofile(pathEnergyIndex)

    return I_total, eIndex
'''
def runTransform2(minrec1 = -np.pi, maxrec1 = np.pi, minrec2 = -np.pi, maxrec2 = np.pi, resolution = sideLength):
    lattice_position = helper.getPositions(fNum)
    lattice_position = np.array(lattice_position).reshape(param[fNum]["atoms"], 3)

    I_total = np.zeros((resolution, resolution, maxEnergyIndex))

    q_x = np.array([0, 0, 0])
    eIndex = np.array([], dtype=np.int32)
    for i in range(resolution):
        q_x = np.array([0, 0, 0])
        q_x = helper.reciprocalPath([minrec1+i/resolution*(maxrec1-minrec1), minrec2, 0], [minrec1+i/resolution*(maxrec1-minrec1), maxrec2,  0], q_x, resolution)
        q_x = np.delete(q_x, 0, axis=0)

        I_aa = getFasterTransform(x[0,:,:],y[0,:,:],z[0,:,:],q_x, 0)
        #I_aa = scatter(q_x[-1, :], lattice_position)
        eIndex = np.append(eIndex, np.argmax(I_aa[:]))
        I_total[i,:, :] = I_aa[:,:maxEnergyIndex]
        print(f'q_x {i} total progress {i/resolution*100:.2f} %', end='\r')
    
    I_total.tofile(pathFourier)
    eIndex.tofile(pathEnergyIndex)

    return I_total, eIndex

if recalculate:
    start = time.time()
    I_total, eIndex = runTransform2(minq1, maxq1, minq2, maxq2, resol)
    print(f'duration: {time.time()-start}')
else:
    I_total = np.fromfile(pathFourier)
    I_total = I_total.reshape((resol, resol, maxEnergyIndex))
    eIndex = np.fromfile(pathEnergyIndex, dtype=np.int32)
print(eIndex)
eIndex = np.unique(eIndex)
eIndex  =np.where(eIndex > np.size(energies),np.size(energies)-1, eIndex)
eIndex = eIndex[energies[eIndex].argsort()]
#values, bins, _ = plt.hist(energies[eIndex], bins=144)
#print(energies[eIndex])


SpecificEnergy = 1
EnergyVariance=0.2
SpecificEnergyIndex =np.argmax(energies > SpecificEnergy)
SpecificEnergyIndexmin =np.argmax(energies > SpecificEnergy-EnergyVariance)
SpecificEnergyIndexmax =np.argmax(energies > SpecificEnergy+EnergyVariance)
extent = [minq1, maxq1, minq2, maxq2]
print(eIndex,energies)
print(SpecificEnergyIndex)
'''
fig, ax = plt.subplots(4, 4)
for i in range(4):
    for j in range(4):
        #I_total[np.where(I_total>1)]=0
        ax[i][j].imshow(np.log(I_total[:, :, int(eIndex[1 * i + j])]), extent=extent)
        ax[i][j].set_aspect(abs(extent[1] - extent[0]) / abs(extent[3] - extent[2]))
        ax[i][j].set_title(f'E = {energies[int(eIndex[1 * i + j])]:.3f}')
# im = plt.imshow(I_total[:,:,eIndex], extent = extent)
plt.tight_layout()
'''
plt.figure()
ttt = np.sum(I_total[:,:,SpecificEnergyIndexmin:SpecificEnergyIndexmax],axis=2)
im = plt.imshow(np.log(ttt), extent = extent)
plt.xticks(xtic,xtlab)
plt.yticks(ytic,ytlab)
plt.ylabel("(-kk0)")
plt.xlabel("(hh0)")
plt.colorbar(label="log $S(Q,\\omega)$")

plt.title(r'T = %i K'%(param[fNum]['temperature']))
# plt.plot(q[1:,0], helper.constants["J_to_meV"] * 4 * param[fNum]["J"]*3.5*(1-np.cos(q[1:,0])), 'w')
#plt.figure()
#im = plt.imshow(I_total[10,:,:], extent = extent)
#maxind = np.where(energies>=180)[0][0]
#plt.imshow(np.log(I_total[10, :,:maxind]).T,aspect = 2/50,origin  ='lower',interpolation='none',vmin =0)
plt.show()
