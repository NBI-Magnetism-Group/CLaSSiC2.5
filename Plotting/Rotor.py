import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()


fig, ax = plt.subplots()
colors = ['r','b','g']
c = 0

fNum = 0
atom = 0
J = param[fNum]["J"]
while fNum < x.shape[0]:
    energy, angle = [], []
    while fNum < x.shape[0] and J == param[fNum]["J"]:
        energy.append(helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"])[np.argmax(np.fft.fft(x[fNum, atom, :]))])
        angle.append(90-np.degrees(np.arccos(z[fNum, atom, -1])))
        fNum = fNum + 1
    
    theoryAngle  = np.linspace(0,90,100)
    theoryEnergy = helper.constants["J_to_meV"] * 4 * param[fNum-1]["J"] * 7/2 * np.sin(np.radians(theoryAngle))

    ax.plot(theoryAngle, theoryEnergy, colors[c], label=f'J = {param[fNum-1]["J"]/helper.constants["boltzmann"]:.0f} K')
    ax.plot(angle, np.abs(energy), colors[c]+'o')
    c = c + 1
    print(np.abs(energy))
    if fNum < x.shape[0]:
        J = param[fNum]["J"]
ax.set_xlabel(r'$\theta$ [degrees]')
ax.set_ylabel(r'Energy [meV]')
ax.legend()
ax.set_title('Rotor')

# time = np.linspace(0,param[fNum][fNum]["dt"]*N,N)
# plt.plot(time,x[:, 0], 'r')
# ax[1].plot(4.135665538536e-12*np.real(np.fft.fftfreq(np.size(x),d=dt)), np.fft.fft(x))
# ax[1].set_xlim(-1,1)
#
# plt.axis('equal')
plt.show()