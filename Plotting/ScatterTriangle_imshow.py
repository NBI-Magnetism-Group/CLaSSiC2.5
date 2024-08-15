import helper
import os, time
import numpy as np
import matplotlib.pyplot as plt
from turtle import shape
from matplotlib import colors
from scipy.signal import find_peaks
from matplotlib import rc
# Enabeling latex
#rc('text', usetex=True)



def getFasterTransform(sx,sy,sz,qScatter, i = 0):
    latticePosition = helper.getPositions(fNum)
    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], i)
    makenew = True #if i want to calculate new fourier

    theta = np.pi/3*0 #rotates data
    sa = np.cos(theta)*sx[:,:]-np.sin(theta)*sy[:,:]
    sb = np.cos(theta)*sy[:,:]+np.sin(theta)*sx[:,:]
    LP2 = np.zeros_like(latticePosition)
    LP2[:,0] = np.cos(theta)*latticePosition[:,0]-np.sin(theta)*latticePosition[:,1]
    LP2[:,1] = np.cos(theta)*latticePosition[:,1]+np.sin(theta)*latticePosition[:,0]
    LP2[:,2] = latticePosition[:,2]

    if not os.path.exists(pathFourier) or makenew:
        print("Calculating")
        start = time.time()
        print("x: ")
        I_total_x = helper.getF_to_I_dp(sa[:,:],qScatter,maxEnergyIndex,param,LP2,0)
        print("y: ")
        I_total_y = helper.getF_to_I_dp(sb[:,:],qScatter,maxEnergyIndex,param,LP2,0)
        print("z: ")
        I_total_z = helper.getF_to_I_dp(sz[:,:],qScatter,maxEnergyIndex,param,LP2,0)
        I_total =   np.sqrt(I_total_x**2 + I_total_y**2 + I_total_z**2)
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
param, x, y, z = helper.getData2(fNum)

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies#/10
maxEnergyIndex = int(param[fNum]["steps"]//2+1)

pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
sideLength = param[fNum]["nUnitCells"]
spin = 2
con_ptr = 1

BField = param[fNum]['magneticField']
Temp = param[fNum]['temperature']
J = param[fNum]['J']*helper.constants['boltzmann']*helper.constants["J_to_meV"]

############## calculate ####################
# The following while loop can be used to delete all the fourier.dat files 

# while os.path.exists(pathFourier):
#     os.remove(pathFourier)
#     fNum += 1
#     pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
# fNum = 0

#Using the helper function to find the q vectors for the triangle

#q_100 = helper.scatterLine(int(sideLength/2))

#q_100 = helper.scatterLine6(int(sideLength/2))
#q_100 = helper.scatterTriangle6(int(sideLength/2))
q_100 = helper.scatterLine5(int(sideLength/2))
nScatter = q_100.shape[0]
q_leng = np.linalg.norm(q_100,axis=1)
# summing over all of the modes (Havn't tested this code with only one simulated mode, but think it should work)

for i in range(x.shape[0]):
    if i ==0:
        I_total = getFasterTransform(x[i,:,:],y[i,:,:],z[i,:,:],q_100, i)
    else:
        I_total += getFasterTransform(x[i,:,:],y[i,:,:],z[i,:,:],q_100, i)


############# Plotting ###############
plt.figure()
maxind = np.where(energies>=400)[0][0] # the max energi on the y axis 

#### find ticks
ytic = np.linspace(0,maxind,9)
ytlab = [str(round(energies[int(ytic[i])],2)) for i in range(len(ytic))]
xtic = np.linspace(0,nScatter,5)
xtlab = [r'$-\pi$',r'$-\pi /2$','0',r'$\pi /2$',r'$\pi$']
#xtlab = [r'$-2\pi$',r'$-\pi$','0',r'$\pi$',r'$2\pi$']

#xtic = np.array([0, (nScatter-1)/(1+3**(3/2)/4), (nScatter-1)*(1+np.sqrt(3)/2)/(1+3**(3/2)/4), nScatter-1])
#xtlab = ['M','$\Gamma$','K','M']

#xtic = np.array([0, (nScatter-1)/3, (nScatter-1)/2, (nScatter-1)*2/3, nScatter-1])
#xtlab = ['$\Gamma$',"$\Gamma$'", "M'","$\Gamma$'",'$\Gamma$']



# aspect: controles the dimensions of the image produced
plt.imshow(np.log(I_total[:,:maxind]).T,aspect = "auto",origin  ='lower',interpolation='none',vmin=0, vmax=12)#,vmax=8,vmin=6)#cmap="nipy_spectral",
plt.yticks(ytic,ytlab)
plt.xticks(xtic,xtlab)
cbar = plt.colorbar()
cbar.set_label(r'Simulated log S($\mathbf{Q}$,$\omega$) [A.U.]')

theoscaller=20
###### Plot theory ######
#calculatets the theoretical value of the q vectors
if param[fNum]['J']>0:
    theo = np.abs(helper.constants["J_to_meV"]*(8*param[fNum]["J"]*spin*(1-np.cos(q_leng[:]))+helper.constants["gFactor"]*helper.constants["bohrMagneton"]*(param[fNum]["magneticField"][-1]+param[fNum]["anisotropyStrength"])))
elif param[fNum]['J']<0: #change theory for AFM
    #a = theoscaller*2*spin*abs(param[fNum]["J"])*np.sqrt(4-(np.cos(q_100[:,0])+1)**2)
    #a = 4*spin*abs(param[fNum]["J"])*abs(np.sin(q_100[:,0]))
    a = theoscaller*2*spin*abs(param[fNum]["J"])*np.sqrt(9-(np.cos(q_100[:,0])+2)**2)
    #a = theoscaller*2*spin*abs(param[fNum]["J"])*np.sqrt((2.1)**2-(np.cos(q_100[:,0])+1.1)**2)
    #gamma = 1/3*(np.cos(q_100[:,0])+2*np.cos(q_100[:,0]/2)*np.cos(q_100[:,1]*np.sqrt(3)/2))
    #a = 3*abs(param[fNum]["J"])*spin*np.sqrt((1-gamma)*(1+2*gamma))
    #a = 10*2*spin*abs(param[fNum]["J"])*np.sqrt(1-(np.cos(q_100[:,0]))**2)
    #a = 10*theoscaller*2*spin*abs(param[fNum]["J"])*np.sqrt((1+0.1)**2-(np.cos(q_100[:,0])+0.1)**2)
    #a = theoscaller*2*spin*abs(param[fNum]["J"])*np.sqrt(1-(np.cos(q_100[:,0])**2))
    """
    #a1 = (2*4*param[fNum]["J"]*spin)**2*(1-np.cos(q_leng[:])**2)
    #a1 = (param[fNum]["J"]*spin)**2*(4**2-2**2*(np.cos(q_leng[:]/np.sqrt(2))*2)**2)*3.5**2
    #a1 = (param[fNum]["J"]*spin)**2*(6**2-2**2*(np.cos(q_100[:,0])+np.cos(q_100[:,0]/2+q_100[:,1]*np.sqrt(3)/2)+np.cos(-q_100[:,0]/2+q_100[:,1]*np.sqrt(3)/2))**2)
    gamma = 1/3*(np.cos(q_100[:,0])+2*np.cos(q_100[:,0]/2)*np.cos(q_100[:,1]*np.sqrt(3)/2))
    a1 = 3.5**2*(3*param[fNum]["J"]*spin)**2*(1-gamma)*(1+2*gamma)
    
    a2 = -8*param[fNum]["J"]*spin*helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["anisotropyStrength"]
    a3 = (helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["anisotropyStrength"])**2
    b = a1+a2+a3
    """
    #theo = helper.constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)+helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1])
    theo = helper.constants["J_to_meV"]*(a+helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1])
    print(theo)
    """
    if param[fNum]["magneticField"][-1] != 0:
        #theo_2 = np.abs(helper.constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)-helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1]))
        theo_2 = np.abs(helper.constants["J_to_meV"]*(a-helper.constants["gFactor"]*helper.constants["bohrMagneton"]*param[fNum]["magneticField"][-1]))
        for i in range(len(theo)):
            if i ==0:
                ind_2 = np.array([np.where(energies>=theo_2[i])[0][0]])
                ind_x = np.array([i])
            else:
                ind_2 = np.append(ind_2,np.where(energies>=theo_2[i])[0][0])
                ind_x = np.append(ind_x,i)
        plt.plot(ind_x,ind_2, '--w',alpha = 0.5)
     
# Puts the theoretical values into arrays, so they can be plotted ontop of the image
"""
for i in range(len(theo)):
    if i ==0:
        ind = np.array([np.where(energies>=theo[i])[0][0]])
        ind_x = np.array([i])
    else:
        ind = np.append(ind,np.where(energies>=theo[i])[0][0])
        ind_x = np.append(ind_x,i)
#print(ind, ind_x)
plt.plot(ind_x,ind*0.05, '--w',alpha = 0.5,label="Theory")
plt.plot(ind_x,ind, '--r',alpha = 0.5,label="Theory times %d"%(theoscaller))
plt.legend()
plt.title(r'T = %i K'%(param[fNum]['temperature']))
#print(energies[ind])
print(param[0])
plt.xlabel(r'$\mathbf{Q}$ = ($\mathit{h00})$ [$Å ^{-1}$]')
plt.ylabel(r'$\mathit{Energy}$ [meV]')
plt.show()

