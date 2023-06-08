import os
import numpy as np
import matplotlib.pyplot as plt
import helper

param, x, y, z = helper.getData()

print('plotting...')
fig, ax = plt.subplots()

colors = ['b','g','r','c','m','y','b']
c = 0
B = param[0]["magneticField"][-1]
i = 0
while i<x.shape[0]:
    zArray = []
    temperatureArray = []
    print(i)
    while i<x.shape[0] and B == param[i]["magneticField"][-1]:
        print(f'B {param[i]["magneticField"][-1]}')
        temperatureArray.append(param[i]["temperature"])
        zArray.append(np.mean(z[i,0,int(param[i]["steps"]/2):]))
        i = i+1

    y = 2.002 * 7/2 * param[i-1]["magneticField"][-1]*9.274009994e-24/1.38064852e-23
    temperature = np.linspace(1,param[-1]["temperature"],50)
    y = y/np.array(temperature)
    langevin = np.cosh(y)/np.sinh(y) - 1/y
    print(f'Magnetic field: {param[i-1]["magneticField"][-1]}')

    plt.plot(temperature, langevin , colors[c], label=f'B = {param[i-1]["magneticField"][-1]}')
    plt.plot(temperatureArray, zArray,colors[c]+'.')
    c = c + 1

    if i<x.shape[0]:
        B = param[i]["magneticField"][-1]

ax.set_xlabel('$T$ [K]')
ax.set_ylabel(r'${\langle s^z \rangle}/{s}$ [-]')
ax.legend()
#
# plt.axis('equal')
plt.show()