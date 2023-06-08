import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0
atom = 0
rows = 3
columns = 4

fig, ax = plt.subplots(rows,columns)
end = int(1e2)
time = np.linspace(0,param[fNum]["dt"]*param[fNum]["steps"],param[fNum]["steps"])
for i in range(rows):
    for j in range(columns):
        if rows > 1:
            ax[i][j].plot(time[:end],x[fNum, 3*i+j, :end], 'r')
            ax[i][j].plot(time[:end],y[fNum, 3*i+j, :end], 'b')
            # ax[i][j].plot(time[:end],z[fNum, 3*i+j, :end], 'g')
        else:
            ax[i][j].plot(time[:end],x[fNum, 3*i+j, :end], 'r')
            ax[i][j].plot(time[:end],y[fNum, 3*i+j, :end], 'b')
            ax[i][j].plot(time[:end],z[fNum, 3*i+j, :end], 'g')
        precession = helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"])[np.argmax(np.fft.fft(x[fNum, atom, :]))]
        print(f'Precession atom {3*i+j}: {precession}')
        print(f'start atom {3*i+j}: {x[fNum, 3*i+j, 0]}, {y[fNum, 3*i+j, 0]}, {z[fNum, 3*i+j, 0]}')

plt.show()
