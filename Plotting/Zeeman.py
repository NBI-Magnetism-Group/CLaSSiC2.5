import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum  = 0
atom = 0

print(x.shape)
fig, ax = plt.subplots(1,3,constrained_layout=True)
time = np.linspace(0,param[fNum]["steps"]*param[fNum]["dt"]*1e12,param[fNum]["steps"])
#energy = data[3::N]


print(f'z difference: {np.max(z[fNum, atom, :])-np.min(z[fNum, atom, :])}')
print(f'final length: {np.sqrt(x[fNum, atom, -1]**2 + y[fNum, atom, -1]**2 + z[fNum, atom, -1]**2)}')
print(f'Average Z-Pos: {np.average(z[fNum, atom, int(len(z[fNum, atom,:]) * 0.9):])}')
precession = helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"])[np.argmax(np.fft.fft(x[fNum, atom, :]))]
theory = helper.constants["Hz_to_meV"]*helper.constants["gamma"]*param[fNum]["magneticField"][-1]/(2*np.pi)
print(f'precession energy: {precession} meV, theoretical: {theory} meV, error: {theory/precession}')
#"""
end = int(2e5)
ax[0].plot(time[:end], x[fNum, atom, :end],'r', label='x')
ax[0].plot(time[:end], y[fNum, atom, :end],'b', label='y')
ax[0].plot(time[:end], z[fNum, atom, :end],'g', label='z')
ax[0].legend(loc='upper right')
#ax[0].set_xlim(0,17.85) #ca. 10 perioder
ax[0].set_xlabel('Time [ps]')
ax[0].set_ylabel('spin [-]')
ax[0].set_title('Spin evolution')

fft_x_limit = 1.5
ax[1].plot(helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"]), np.fft.fft(x[fNum, atom, :]).real/np.max(np.fft.fft(x[fNum, atom, :]).real))
ax[1].set_xlim(-fft_x_limit, fft_x_limit)
ax[1].set_xlabel('Energy [meV]')
ax[1].set_ylabel('relative counts [-]')
ax[1].set_title('Fourier transform')

ax[2].plot(x[fNum, atom, :end], y[fNum, atom, :end],'r')
ax[2].set_xlabel('x')
ax[2].set_ylabel('y')
ax[2].set_title('x-y plane')

# plt.axis('equal')
plt.show()
#"""