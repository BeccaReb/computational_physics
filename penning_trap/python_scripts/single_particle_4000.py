import numpy as np
import matplotlib.pyplot as plt

Penningtrap_data_1 =np.genfromtxt("../textfiles/single_particle_1_4000_rk.txt", dtype="float", delimiter=",", skip_header=0)


time = Penningtrap_data_1[:,0]
z = Penningtrap_data_1[:,3]


#--------------------------------------------------------------------------------------
#plots particle 1
fig, ax = plt.subplots(figsize = (8, 6))
ax.plot(time, z, color = "#ee5d5f", label = "Particle 1")
ax.scatter(time[0], z[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
#ax.set_title("Single particle z-position as a function of time" , fontsize = 16)
ax.set_xlabel("time [$\mu$s]", fontsize = 18)
ax.set_ylabel("z [$\mu$m]", fontsize = 18)
ax.set_xlim(-2, 59)
ax.set_ylim(-23, 22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(color = "#e7e7e7", linestyle = "dashed")
plt.legend(fontsize = 13, loc = "upper right")
#plt.savefig("../figures/single_particle_z.pdf")
plt.show() 
#--------------------------------------------------------------------------------------

