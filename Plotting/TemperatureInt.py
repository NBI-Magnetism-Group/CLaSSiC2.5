import os, time
import numpy as np
import matplotlib.pyplot as plt
import helper
from scipy import optimize as op

param, x, y, z = helper.getData2()
sideLength = param[0]["nUnitCells"]
size = int(sideLength)

frequencies = np.fft.fftfreq(param[0]["steps"], param[0]["dt"])
energies = 4.1357e-12 * frequencies/np.pi**2
EnergyVariance = 0.2
qvariance = 0.05
SpecificEnergyM = 5
SpecificEnergyIndexMmax = np.argmax(energies > SpecificEnergyM+EnergyVariance)
SpecificEnergyIndexMmin = np.argmax(energies > SpecificEnergyM-EnergyVariance)
SpecificEnergyG = 1
SpecificEnergyIndexGmax = np.argmax(energies > SpecificEnergyG+EnergyVariance)
SpecificEnergyIndexGmin = np.argmax(energies > SpecificEnergyG-EnergyVariance)

print(SpecificEnergyIndexMmax,energies[SpecificEnergyIndexMmax])
print(SpecificEnergyIndexMmin,energies[SpecificEnergyIndexMmin])
print(SpecificEnergyIndexGmax,energies[SpecificEnergyIndexGmax])
print(SpecificEnergyIndexGmin,energies[SpecificEnergyIndexGmin])

temperatureArray = []
scalling1 = np.pi/3*8
scalling2 = np.pi/3*8*np.sqrt(3)/2
q_G = np.array([0, 0, 0])
q_M = np.array([0, 0, 0])
qminM = -0.25
qmaxM = 0.25
qminG = -0.25
qmaxG = 0.25

deltaM = np.linspace(qminM,qmaxM,size)
deltaG = np.linspace(qminG,qmaxG,size)
deltaMH = np.linspace(qminM,qmaxM,100)
deltaGH = np.linspace(qminG,qmaxG,100)

#q_G = helper.reciprocalPath([scallingy/4, scallingx/(np.sqrt(3)*4), 0], [scallingy/4, scallingx*(1-1/(np.sqrt(3)*4)), 0], q_G, size)
#q_G = helper.reciprocalPath([scallingy/4, scallingx*qmin, 0], [scallingy/4, scallingx*qmax, 0], q_G, size)
#q_G = helper.reciprocalPath([scallingy/4, 0, 0], [scallingy*3/4, 0, 0], q_G, size)
#q_G = helper.reciprocalPath([scallingx*(1.5+qmin)/2, scallingx*(-2*qmin), 0], [scallingx*(1.5+qmax)/2, scallingx*(-2*qmax), 0], q_G, size)
q_G = helper.reciprocalPath([scalling1*(1+qminG)/2, scalling2*(-1+qminG), 0], [scalling1*(1+qmaxG)/2, scalling2*(-1+qmaxG), 0], q_G, size)

#q_M = helper.reciprocalPath([0, scallingx/(np.sqrt(3)*4), 0], [0, scallingx*(1-1/(np.sqrt(3)*4)), 0], q_M, size)
q_M = helper.reciprocalPath([scalling1*(1.5+qminM)/2, scalling2*(-2*qminM), 0], [scalling1*(1.5+qmaxM)/2, scalling2*(-2*qmaxM), 0], q_M, size)

q_G = np.delete(q_G, 0, axis=0)
q_M = np.delete(q_M, 0, axis=0)
print('plotting...')

colors = ['b','g','r','c','m','y','b']
c = 0
B = param[0]["magneticField"][-1]
i = 0
maxEnergyIndex = int(param[0]["steps"]//2+1)

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

def Lorentzian(x, amp, cen, wid):
    return (amp*wid**2/((x-cen)**2+wid**2))

M_sums =np.empty([size,1])
G_sums =np.empty([size,1])
M_amp = np.array([])
M_cen = np.array([])
M_wid = np.array([])
G_amp = np.array([])
G_cen = np.array([])
G_wid = np.array([])
for i in range(helper.countData()):
    param, x, y, z = helper.getData2(i)
    I_totalM = getFasterTransform(x[0,:,:],y[0,:,:],z[0,:,:],q_M, i)
    I_totalG = getFasterTransform(x[0,:,:],y[0,:,:],z[0,:,:],q_G, i)

    M_sum=np.sum(I_totalM[:,SpecificEnergyIndexMmin:SpecificEnergyIndexMmax],axis=1)
    G_sum=np.sum(I_totalG[:,SpecificEnergyIndexGmin:SpecificEnergyIndexGmax],axis=1)


    M_sums = np.append(M_sums,np.array([M_sum]).T,axis=1)
    G_sums = np.append(G_sums,np.array([G_sum]).T,axis=1)

    M_val, M_cov = op.curve_fit(Lorentzian,deltaM,M_sum,p0=(np.max(M_sum),0,qmaxM),bounds=((0,qminM,0),(np.max(M_sum)*10,qmaxM,2)))
    M_amp = np.append(M_amp,M_val[0])
    M_cen = np.append(M_cen,M_val[1])
    M_wid = np.append(M_wid,M_val[2])

    G_val, M_cov = op.curve_fit(Lorentzian,deltaG,G_sum,p0=(np.max(G_sum),0,qmaxG),bounds=((0,qminG,0),(np.max(G_sum)*10,qmaxG,2)))
    G_amp = np.append(G_amp,G_val[0])
    G_cen = np.append(G_cen,G_val[1])
    G_wid = np.append(G_wid,G_val[2])

    temperatureArray.append(param[0]["temperature"])
    



fig1, ax1 = plt.subplots(5, 4)
for i in range(5):
    for j in range(4):
        ax1[i,j].plot(deltaM,M_sums[:,1+i*4+j],'o',label="%i K"%(temperatureArray[i*4+j]))
        ax1[i,j].plot(deltaMH,Lorentzian(deltaMH,M_amp[i*4+j],M_cen[i*4+j],M_wid[i*4+j]))
        ax1[i,j].legend()
        if j == 0:
            ax1[i,j].set_ylabel("Intetnsity [a.u.]")
        if i==4:
            ax1[i,j].set_xlabel("q=(-2$\Delta$, 1.5+$\Delta$,0)")
 


fig2, ax2 = plt.subplots(5, 4)
for i in range(5):
    for j in range(4):
        ax2[i,j].plot(deltaG,G_sums[:,1+i*4+j],'o',label="%i K"%(temperatureArray[i*4+j]))
        ax2[i,j].plot(deltaGH,Lorentzian(deltaGH,G_amp[i*4+j],G_cen[i*4+j],G_wid[i*4+j]))
        ax2[i,j].legend()
        if j == 0:
            ax2[i,j].set_ylabel("Intetnsity [a.u.]")
        if i == 4:
            ax2[i,j].set_xlabel("q=(-1+$\Delta$, 1+$\Delta$,0)")

fig3, ax3 = plt.subplots()
ax3.set_xlabel("T [K]")
ax3.set_ylabel("S($\Gamma$', 1 meV) [a.u.]", color="C0")
p0 = ax3.plot(temperatureArray,G_amp,'o', color="C0", label="$\Gamma$'")
ax3.tick_params(axis='y',labelcolor="C0")

ax4 =ax3.twinx()
ax4.set_ylabel("S(M', 5 meV) [a.u.]", color="C1")
p1 = ax4.plot(temperatureArray,M_amp,'o', color="C1", label="M'")
ax3.axvline(71,color="grey")
ax3.text(80,np.max(G_amp)*0.8,'T$_N$', color="grey")
ax4.tick_params(axis='y',labelcolor="C1")
axs = p0+p1
labs = [l.get_label() for l in axs]
ax3.legend(axs,labs,loc=0)

plt.figure()
plt.plot(temperatureArray,1/G_wid,'o',label="$\Gamma$'")
plt.plot(temperatureArray,1/M_wid,'o',label="M'")
plt.axvline(71,color="grey")
plt.text(80,np.max(1/G_wid)*0.8,'T$_N$', color="grey")
plt.xlabel("T [K]")
plt.ylabel("$\\xi $ [Å]")
plt.legend()

plt.show()