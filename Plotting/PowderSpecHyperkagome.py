import helper
import os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def construct_q(size, max):
    q = np.array([0, 0, 0])
    a = 2*np.pi/(size*8)
    for k in range(size*max+1):
        for j in range(size*max+1):
            for i in range(size*max+1):
                q = np.vstack((q, np.array([a * ( i ), a * ( j ), a * ( k )])))

    q = np.delete(q, 0, axis=0)
    return q

def construct_max_q_space(size, max, max_vec):
    # Gives all q vectores shorter then that max_vec 
    max_vec = np.array(max_vec)
    q = np.array([0, 0, 0])
    a = 2*np.pi/(size*8)
    for k in range(size*max+1):
        for j in range(size*max+1):
            for i in range(size*max+1):
                q = np.vstack((q, np.array([a * ( i ), a * ( j ), a * ( k )])))

    q = np.delete(q, 0, axis=0)
    max_leng = np.linalg.norm(a*size * max_vec)
    print(f"max_len: {max_leng}")
    len_arr = np.linalg.norm(q, axis=1)
    print(f"Deleted q len: {np.min(len_arr[np.where(len_arr>max_leng)[0]])}")
    q = np.delete(q,np.where(len_arr>max_leng)[0],axis=0)
    len_arr_new = np.delete(len_arr,np.where(len_arr>max_leng)[0],axis=0)
    print(f"max q: {q[np.argmax(len_arr_new)]}")
    return q, len_arr_new

def construct_q_1D(size, direction, max):
    q = [0,0,0]
    a = 2*np.pi/(size*8)
    for i in range(1,size*max+1):
        q_new = np.array([direction[0]*a*i,direction[1]*a*i,direction[2]*a*i])
        q = np.vstack((q, q_new))
    return q

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

def powder_spect(x,y,z,size, max, max_vec,nBins,Filenum):
    # Calculates the powder averaged intensities
    q, leng = construct_max_q_space(size, max, max_vec)
    I = getFasterTransform(x,y,z,q,Filenum)
    I_M = np.zeros((nBins,I.shape[1]))
    bin_edg = np.histogram_bin_edges(leng, bins = nBins-1)
    print(f'bin_edg: {bin_edg}')
    Bin_ind = np.digitize(leng,bin_edg)
    for i in range(int(np.max(Bin_ind)-1)):
        I_M[i,:] = np.mean(I[np.where(Bin_ind==(i+1)),:],axis=1)
    
    return bin_edg, I_M


param, x, y, z = helper.getData()
fNum = 0
sideLength = param[fNum]["nUnitCells"]
latticePosition = helper.getPositions(fNum)
scatterLength = int(sideLength)
frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = int(param[fNum]["steps"]//2+1)

#This is the conversion factor to real units, in the simulation the cubic sidelength is 8, while in reality it's 12.39 Å
con_ptr = 8/12.39

############## fourier10 = powder spec (reminder for the filename) #################
# Be carefull, takes a long time!!!
# Gives the powder averaged intensities and the bin edges 
bin_edg, I_M = powder_spect(x,y,z,scatterLength,3,[3,0,0],35,10)

############## Plotting #################
plt.figure()
maxind = np.where(energies>=0.6)[0][0]
ytic = np.linspace(0,maxind,7)
ytlab = [str(round(energies[int(ytic[i])],5)) for i in range(len(ytic))]
xtic = np.arange(1,len(bin_edg),2)
xtlab = [str(round(np.linalg.norm(bin_edg[int(j)])*con_ptr,2)) for j in xtic]

plt.imshow(np.log(I_M[:,:maxind]).T,aspect = 1/500,origin  ='lower',interpolation='none',cmap="gist_ncar",vmin=8,vmax = 20)#,vmax=12)#,vmax=8,vmin=6)#cmap="nipy_spectral",
plt.yticks(ytic,ytlab)
plt.xticks(xtic,xtlab)
plt.colorbar()
plt.xlabel('Q (Å^-1)')
plt.ylabel('E (meV)')
plt.show()
