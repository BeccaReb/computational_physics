import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

data1 = np.genfromtxt("../textfiles/particle_1_4000_with_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)
data2 = np.genfromtxt("../textfiles/particle_2_4000_with_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)

x1 = data1[:,1]
y1 = data1[:,2]
z1 = data1[:,3]

x2 = data2[:,1]
y2 = data2[:,2]
z2 = data2[:,3]

fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.grid()

ax.scatter(x1[0], y1[0], z1[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
ax.scatter(x2[0], y2[0], z2[0], zorder = 3, color = "#3d8c40", label = "Start particle 2")
ax.plot(x1, y1, z1, c = "#ee5d5f", label = "Particle 1")
ax.plot(x2, y2, z2, c = "#01377d", label = "Particle 2")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

font = {'size': 20}
ax.tick_params('z', labelsize=font['size'])
plt.legend(fontsize = 18)
#ax.set_title('3D Scatter Plot')
#ax.set_zticks([])

plt.legend(fontsize = 18)

# Set axes label
ax.set_xlabel('x [$\mu$m]', labelpad=22, fontsize = 18)
ax.set_ylabel('y [$\mu$m]', labelpad=22, fontsize = 18)
ax.set_zlabel('z [$\mu$m]', labelpad=22, fontsize = 18)

#plt.savefig("../figures/3d_plot_two_particles_with_interaction.pdf")
plt.show()