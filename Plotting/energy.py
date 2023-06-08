import helper
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

param, x, y, z = helper.getData()
fNum = 0

fig, ax = plt.subplots()
end = int(1e6)
time = np.linspace(0,param[fNum]["dt"]*param[fNum]["steps"],param[fNum]["steps"])
energy = helper.getEnergy(fNum)
ax.plot(time[1:end], energy[1:end]/param[fNum]["atoms"])
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
print(np.isnan(np.sum(energy)))
print(energy[0])
# ax.set_yscale("log")
ax.set_xlabel("time [s]")
ax.set_ylabel("energy/N [meV]")

plt.show()
