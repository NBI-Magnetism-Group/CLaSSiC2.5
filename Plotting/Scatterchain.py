import helper
import os, time
import numpy as np
import matplotlib.pyplot as plt
from turtle import shape
from matplotlib import colors
from scipy.signal import find_peaks
from matplotlib import rc
# Enabeling latex
rc('text', usetex=True)



def getFasterTransform(sx,sy,sz,qScatter, i = 0):
    latticePosition = helper.getPositions(fNum)
    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], i)
    if not os.path.exists(pathFourier):
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

def prepdata(I_plot,used_q,usePeaks = True):
    nScatter = len(used_q)
    xData, yData, zData = [], [], []
    maxF = np.log10(np.max(I_plot[1:]))
    for i in range(nScatter):
            if usePeaks:
                peaks, _ = find_peaks(I_plot[i, :], height=1e5, distance=2000)
                # Find peaks cannot check if the peaks at index 0
                if I_plot[i, 0] >= np.max(I_plot[i, :])/10:
                    yData.append(np.abs(energies[0]))
                    xData.append(used_q[i])
                    zData.append(max([np.log10(np.max(I_plot[i,0]))-maxF, -5]))
                # If there are no peaks take the highest value
                if peaks.size == 0:
                    yData.append(np.abs(energies[np.argmax(I_plot[i,:])]))
                    xData.append(used_q[i])
                    zData.append(max([np.log10(np.max(I_plot[i,:]))-maxF, -5]))
                    # Plot all the peaks
                else:
                    for peak in peaks:
                        yData.append(np.abs(energies[peak]))
                        xData.append(used_q[i])
                        zData.append(max([np.log10(np.max(I_plot[i, peak]))-maxF, -5]))
            else:
                yData.append(np.abs(energies[np.argmax(I_plot[i,:])]))
                xData.append(used_q[i])
                zData.append(max([np.log10(np.max(I_plot[i,:]))-maxF, -5]))        
    return xData, yData, zData

def construct_q_1D(size, direction, max):
    q = [0,0,0]
    a = 2*np.pi/(size*8)
    for i in range(1,size*max+1):
        q_new = np.array([direction[0]*a*i,direction[1]*a*i,direction[2]*a*i])
        q = np.vstack((q, q_new))
    
    return q

def get_avg_M(q,I,maxEnergyIndex,windowith=100,plot_int=False,scale=False):
    for i in range(q.shape[0]):
        if plot_int:
            plt.figure()
            helper.plotSmooth(plt, energies[:maxEnergyIndex-1], I[i,:maxEnergyIndex])
            plt.yscale('log')
            plt.title(str(q[i]*con_ptr))
            plt.xlabel('E[meV]')
            plt.ylabel('I')
        xSm, ySm = helper.getSmooth(energies[:maxEnergyIndex-1], I[i,:maxEnergyIndex],windowith)
        if scale:
            ySm = ySm/np.max(ySm)
        if i==0:
            ySmt = ySm
        else:
            ySmt = np.vstack((ySmt,ySm))

    return xSm,ySmt


############# import and set variables ##############
fNum = 0
param, x, y, z = helper.getData()

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = int(param[fNum]["steps"]//2+1)

pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
sideLength = param[fNum]["nUnitCells"]

con_ptr = 1

BField = param[fNum]['magneticField']
Temp = param[fNum]['temperature']
J = param[fNum]['J']*helper.constants['boltzmann']*helper.constants["J_to_meV"]

q_100 = helper.scatterLine(int(sideLength/2))
nScatter = q_100.shape[0]
q_leng = np.linalg.norm(q_100,axis=1)

for i in range(x.shape[0]):
    if i == 0:
        I_total = getFasterTransform(x[i,:,:],y[i,:,:],z[i,:,:],q_100, i)
    else:
        I_total += getFasterTransform(x[i,:,:],y[i,:,:],z[i,:,:],q_100, i)

I_total = I_total.real/sideLength # type : ignore

############# Plotting ###############
plt.figure()
maxind = np.where(energies>=4)[0][0]

#### find ticks
ytic = np.linspace(0,maxind,9)
ytlab = [str(round(energies[int(ytic[i])],2)) for i in range(len(ytic))]
xtic = np.array([0,5,10,15,20])
xtlab = [r'$-\pi$',r'$-\pi /2$','0',r'$\pi /2$',r'$\pi$']

plt.imshow(np.log(I_total[:,:maxind]).T,aspect = 1/500,origin  ='lower',interpolation='none',vmin =0,vmax=7)#,vmax=12)#,vmax=8,vmin=6)#cmap="nipy_spectral",
plt.yticks(ytic,ytlab)
plt.xticks(xtic,xtlab)
cbar = plt.colorbar()
cbar.set_label(r'Intensity [A.U.]')

###### Plot theory ######
if param[fNum]['J']>0:
    theo = np.abs(helper.constants["J_to_meV"]*(4*param[fNum]["J"]*3.5*(1-np.cos(q_leng[:]))+helper.constants["gFactor"]*helper.constants["bohrMagneton"]*(param[fNum]["magneticField"][-1]+param[fNum]["anisotropyStrength"])))
else:
    a1 = (4*param[fNum]["J"]*3.5)**2*(1-np.cos(q_leng[:])**2)
    a2 = -8*param[fNum]["J"]*3.5*helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["anisotropyStrength"]
    a3 = (helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["anisotropyStrength"])**2
    b = a1+a2+a3
    theo = helper.constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)+helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1])
    if param[fNum]["magneticField"][-1] != 0:
        theo_2 = np.abs(helper.constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)-helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1]))
        for i in range(len(theo)):
            if i ==0:
                ind_2 = np.array([np.where(energies>=theo_2[i])[0][0]])
                ind_x = np.array([i])
            else:
                ind_2 = np.append(ind_2,np.where(energies>=theo_2[i])[0][0])
                ind_x = np.append(ind_x,i)
        plt.plot(ind_x,ind_2, '--w',alpha = 0.5)
     
for i in range(len(theo)):
    if i ==0:
        ind = np.array([np.where(energies>=theo[i])[0][0]])
        ind_x = np.array([i])

    else:
        ind = np.append(ind,np.where(energies>=theo[i])[0][0])
        ind_x = np.append(ind_x,i)

print(ind, ind_x)
plt.plot(ind_x,ind, '--w',alpha = 0.5)
plt.title(r'T = 0K, B = 3T, D = 0.1T')
print(energies[ind])
print(param[0])
plt.xlabel(r'\textbf{Q} = (\textit{h00}) [$Å ^{-1}$]')
plt.ylabel(r'\textit{Energy} [meV]')
plt.show()