import helper
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

param, x, y, z = helper.getData()
fNum  = 0
atom = 0
nAtoms = param[fNum]["atoms"]
distance = 0.001

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# line, = ax.plot([atom*distance, atom*distance+x[fNum, atom, 0]], [0, y[fNum, atom, 0]],[0,z[fNum, atom, 0]/100],lw=2,color="b",marker='o')
lines = []

for i in range(nAtoms):
    lobj = ax.plot([atom*distance, atom*distance+x[fNum, atom, 0]], [0, y[fNum, atom, 0]],[0,z[fNum, atom, 0]/1000],lw=2)[0]
    lines.append(lobj)

def animate(i):
    ax.set_title(f"time: {i} ps, length: {x[fNum, 0, i]**2 + y[fNum, 0, i]**2:.2f}")
    for atom, line in enumerate(lines):
        line.set_data([np.array([atom*distance, atom*distance+x[fNum, atom, i]]), np.array([0,y[fNum, atom, i]])])
        line.set_3d_properties(np.array([0,z[fNum, atom, i]/1000]))
    return lines

ax.set(xlim=([-0.002, distance*(nAtoms)+0.002]),ylim=([-distance*(nAtoms)/2-0.002, distance*(nAtoms)/2+0.002]), zlim=([0,0.001]))
"""Moved together to ax.set()
ax.set_xlim3d([-0.002, distance*(nAtoms)+0.002])
ax.set_ylim3d([-distance*(nAtoms)/2-0.002, distance*(nAtoms)/2+0.002])
ax.set_zlim3d([0,0.001])
"""
anim = animation.FuncAnimation(fig, animate,
                               frames=range(0, x.shape[2], 4), interval=100, blit=False)
plt.show()
