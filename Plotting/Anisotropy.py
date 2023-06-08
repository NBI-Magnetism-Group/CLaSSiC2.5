import numpy as np
import matplotlib.pyplot as plt
import helper


angleArray = []
energyArray = []
param, x, y, z = helper.getData()

for i in range(x.shape[0]):
    angleArray.append(np.degrees(np.arccos(z[i,0,0])))
    energyArray.append(6.241509e21*1.0545718e-34*2*np.pi*np.abs(np.fft.fftfreq(np.size(x[i,0,:]),d=param[i]["dt"])[np.argmax(np.abs(np.fft.fft(x[i,0,:])))]))
    print(angleArray[-1])

fig, ax = plt.subplots()

time = np.linspace(0,np.size(x)*param[0]["dt"],np.size(x))
#energy = data[3::N]

print(f'shape {x.shape}')
print(f'size: {np.size(x)}')
print(f'z difference: {np.max(z)-np.min(z)}')
print(np.size(np.fft.fft(x)))
print('plotting...')

print(f'precession energy: {4.135665538536e-12*np.abs(np.fft.fftfreq(np.size(x),d=param[0]["dt"])[np.argmax(np.abs(np.fft.fft(x[0,0,:])))])} MeV')

ax.plot(angleArray, 6.241509e21*2.002* 9.274009994e-24*10*np.cos(np.radians(angleArray)),'r')
ax.plot(angleArray, energyArray, 'bo')
ax.set_xlabel(r'$\theta$ [degrees]')
ax.set_ylabel(r'Energy [meV]')
ax.set_title("Anisotropy")

# plt.axis('equal')
plt.show()