import numpy as np
import matplotlib.pyplot as plt 
import helper

param, x, y, z = helper.getData()
position = helper.getPositions()
fNum = 0

plt.plot(position[:, 0], position[:, 1], 'o', ms=2.5)
multiplier = 0.5 / np.max(x[fNum, :, 0])
print(f'Multiplier: {multiplier}')
for i in range(position.shape[0]):
    plt.plot([position[i, 0], position[i, 0]+multiplier * x[fNum, i, 0]], [position[i, 1], position[i, 1] + multiplier * y[fNum, i, 0]], 'r')

plt.show()