import helper
import os, time
import numpy as np
import matplotlib.pyplot as plt
from turtle import shape
from matplotlib import colors
from scipy.signal import find_peaks
from matplotlib import rc
#enable latex
rc('text', usetex=True)


def getFasterTransform(x,y,z,qScatter, i = 0):
    #Has helper.py calculate the Fourier trans. and makes them into intensities, an saves them and the i'th fourier file.
    latticePosition = helper.getPositions(fNum)
    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], i)
    if not os.path.exists(pathFourier):
        print("Calculating")
        start = time.time()
        print("x: ")
        I_total_x = helper.getF_to_I_dp(x[fNum,:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        print("y: ")
        I_total_y = helper.getF_to_I_dp(y[fNum,:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        print("z: ")
        I_total_z = helper.getF_to_I_dp(z[fNum,:,:],qScatter,maxEnergyIndex,param,latticePosition,0)
        I_total = I_total_x**2 + I_total_y**2 + I_total_z**2
        I_total.tofile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
        print(f'total duration: {time.time()-start}')
    else:
        print("Pulling fourier data from file")
        I_total = np.fromfile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
    return I_total

def prepdata(I_plot,used_q,usePeaks = True):
    # Preps the data for plotting
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
    # Constructs all q vectors in a certain directrion out to a max vector, times a
    q = [0,0,0]
    a = 2*np.pi/(size*8)
    for i in range(1,size*max+1):
        q_new = np.array([direction[0]*a*i,direction[1]*a*i,direction[2]*a*i])
        q = np.vstack((q, q_new))
    
    return q

def get_avg_M(q,I,maxEnergyIndex,windowith=100,plot_int=False):
    # Used to get the smoothenend intensities
    for i in range(q.shape[0]):
        if plot_int:
            plt.figure()
            helper.plotSmooth(plt, energies[:maxEnergyIndex-1], I[i,:maxEnergyIndex])
            plt.yscale('log')
            plt.title(str(q[i]*con_ptr))
            plt.xlabel('E[meV]')
            plt.ylabel('I')
        xSm, ySm = helper.getSmooth(energies[:maxEnergyIndex-1], I[i,:maxEnergyIndex],windowith)
        if i==0:
            ySmt = ySm
        else:
            ySmt = np.vstack((ySmt,ySm))

    return xSm,ySmt

##### These aren't used any more but can be usefull to investergate the q-vectors 
def construct_q(size, max):
    q = np.array([0, 0, 0])
    a = 2*np.pi/(size*8)
    for k in range(size*max+1):
        for j in range(size*max+1):
            for i in range(size*max+1):
                q = np.vstack((q, np.array([a * ( i ), a * ( j ), a * ( k )])))
    q = np.delete(q, 0, axis=0)
    
    return q

def make_unique(q,I):
    leng = np.linalg.norm(q,axis=1)
    unq = np.unique(leng)
    I_avg = I[0]
    for i in range(len(unq)-1):
        ind = np.where(unq[i+1] == leng)[0]
        I_token = np.mean(I[ind],axis=0)
        I_avg = np.vstack((I_avg,I_token))
    
    return unq,I_avg

def make_unique_ind(q):
    leng = np.linalg.norm(q,axis=1)
    unq = np.unique(leng)
    ind = np.zeros(20)
    for i in range(len(unq)-1):
        token1 = np.zeros(20)
        token2 = np.where(unq[i+1] == leng)[0]
        for i in range(len(token2)):
            token1[i] = token2[i]
        ind = np.vstack((ind,token1))
    return unq,ind


############# import and set variables ##############
fNum = 0
param, x, y, z = helper.getData()

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = int(param[fNum]["steps"]//2+1)

pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
sideLength = param[fNum]["nUnitCells"]

#This is the conversion factor to real units, in the simulation the cubic sidelength is 8, while in reality it's 12.39 Å
con_ptr = 8/12.39

BField = param[fNum]['magneticField']
Temp = param[fNum]['temperature']
J = param[fNum]['J']*helper.constants['boltzmann']*helper.constants["J_to_meV"]

############## calculate ####################
# calculate the diretrion specific q vectors: 
q_100 = construct_q_1D(sideLength,[1,0,0],3)
q_110 = construct_q_1D(sideLength,[1,1,0],3)
q_111 = construct_q_1D(sideLength,[1,1,1],3)
q_arr = np.array([q_100,q_110,q_111])

# and there corrosponding intensities.
I_100 = getFasterTransform(x,y,z,q_100, i = 1)
I_110 = getFasterTransform(x,y,z,q_110, i = 2)
I_111 = getFasterTransform(x,y,z,q_111, i = 3)
I_arr = [I_100,I_110,I_111]

# Just a reminder for the naming of the file 
#fourier1 = (100)
#fourier2 = (110)
#fourier3 = (111)
#fourier10 = powder spect max [200]

titles = ['100','110','111']

# add_to_title = f'B = { param[0]['magneticField'] }, T = {param[0]['temperature']}'

maxEnergyplot = int(np.where(energies>=1.2)[0][0])
maxEnergyplot_im = int(np.where(energies>=0.4)[0][0])
nScatter = q_100.shape[0] # type : ignore

avg = np.array([5]) #The size of the window of intensities to average, may only be whole numbers

asp = 1/200 #aspect of the imshow plot

############## Plotting ####################
for i in range(q_arr.shape[0]):
    xSm, ySmt = get_avg_M(q_arr[i],I_arr[i],maxEnergyplot,avg[0],False) #if you set this to true, you get all the intensity plots
    ###### find ticks ##########
    maxind = np.where(xSm>0.6)[0][0]
    ytic = np.linspace(0,maxind,7)
    ytlab = [str(round(energies[int(ytic[i])],2)) for i in range(len(ytic))]
    xtic = np.arange(1,nScatter,2)
    xtlab = [str(round(np.linalg.norm(q_arr[i,int(j)])*con_ptr,2)) for j in xtic]
    plt.figure()
    if i==0:
        plt.title(r'Scatterplot, $\log_{10}$ scale, [\textit{h00}],T = 0.00935 K, B = 4 T, J = -0.186 K')
    if i==1:
        plt.title(r'Scatterplot, $\log_{10}$ scale, [\textit{hh0}],T = 0.00935 K, B = 4 T, J = -0.186 K')
    if i==2:
        plt.title(r'Scatterplot, $\log_{10}$ scale, [\textit{hhh}],T = 0.00935 K, B = 4 T, J = -0.186 K')
    plt.imshow(np.log(ySmt[:,:maxind]).T,aspect = asp,origin  ='lower',interpolation='none',cmap="gist_ncar",vmax=23,vmin =7)#,vmax=12)#,vmax=8,vmin=6)#cmap="nipy_spectral",
    plt.yticks(ytic,ytlab)
    plt.xticks(xtic,xtlab)
    plt.colorbar()
    plt.xlabel(r'\textbf{Q} [$Å^{-1}$]')
    plt.ylabel(r'\textit{E} [meV]')


# if Do_you_want_imshow:
#     fig, ax = plt.subplots(3,2)
# else:
#     fig, ax = plt.subplots(3)

# for i in range(3):
#     fi, axi = plt.subplots(6,3)
#     q_plot = np.linalg.norm(q_arr[i],axis =1)*(8/12.39)
#     # I_plot = I_arr[i]
#     xSm, I_plot = get_avg_M(q_arr[i],I_arr[i],maxEnergyplot,avg[0],False)

#     xDat, yDat, zDat = prepdata(I_plot,q_plot)

#     if Do_you_want_imshow:
#         ticlab = []
#         tic = []
#         ax[i,0].imshow( np.log(I_plot[:,:maxEnergyplot_im]).T, aspect = 1/25, origin  ='lower', cmap = plt.get_cmap('nipy_spectral'))

#         for j in range(len(q_plot)):
#             if j%2 ==0:

#                 ticlab.append(str(np.round(q_plot[j],2)))
#                 tic.append(j)

#         if i==0:
#             print((tic,ticlab))
        
#         ax[i,0].set_xticks(tic,ticlab)

#         ax[i,1].plot(xDat,yDat, helper.constants["colors"][fNum]+'o')
#         helper.plotTheory(param[fNum]["geometry"], sideLength, param[fNum], ax[i,1], fNum)
#         im = ax[i,1].scatter(xDat,yDat,c=zDat,cmap="Wistia", zorder=10)
#         fig.colorbar(im, ax=ax[i,1])
#         ax[i,1].set_title(titles[i])
    
#     else:
#         ax[i].plot(xDat,yDat, helper.constants["colors"][fNum]+'o')
#         helper.plotTheory(param[fNum]["geometry"], sideLength, param[fNum], ax[i], fNum)
#         im = ax[i].scatter(xDat,yDat,c=zDat,cmap="Wistia", zorder=10)
#         fig.colorbar(im, ax=ax[i])
#         ax[i].set_title(titles[i])

#     print(xSm.shape,I_plot.shape)
#     token = 1
#     for j in range(3):
#         for k in range(6):


#             helper.plotPeaks(axi[k,j], xSm, I_plot[token,:])
#             axi[k,j].plot(xSm, I_plot[token,:], c= helper.constants["colors"][fNum])
#             # helper.plotSmooth(axi[k,j], xSm, I_plot[token,:])
#             axi[k,j].set_yscale("log")
#             # axi[k,j].set_title(f"q = {np.round(q_plot,4)}")

#             if token < (nScatter-1):
#                 token +=1

plt.show()