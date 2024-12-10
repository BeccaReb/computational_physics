
import numpy as np
import matplotlib.pyplot as plt

#no particle interaction case
Penningtrap_data_1_no_interaction =np.genfromtxt("../textfiles/particle_1_4000_no_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)
Penningtrap_data_2_no_interaction =np.genfromtxt("../textfiles/particle_2_4000_no_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)

#particle interaction case
Penningtrap_data_1_interaction =np.genfromtxt("../textfiles/particle_1_4000_with_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)
Penningtrap_data_2_interaction =np.genfromtxt("../textfiles/particle_2_4000_with_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)

time = Penningtrap_data_1_no_interaction[:,0]

#no particle interaction case
x_particle_1_no_interaction = Penningtrap_data_1_no_interaction[:,1]
y_particle_1_no_interaction = Penningtrap_data_1_no_interaction[:,2]
x_particle_2_no_interaction = Penningtrap_data_2_no_interaction[:,1]
y_particle_2_no_interaction = Penningtrap_data_2_no_interaction[:,2] 

#particle interaction case
x_particle_1_interaction = Penningtrap_data_1_interaction[:,1]
y_particle_1_interaction = Penningtrap_data_1_interaction[:,2]
x_particle_2_interaction = Penningtrap_data_2_interaction[:,1]
y_particle_2_interaction = Penningtrap_data_2_interaction[:,2]

##-------------------------------------------------------------------
#plots the two particles with no interactions 
fig, ax = plt.subplots(figsize = (8,8))
ax.scatter(x_particle_1_no_interaction[0], y_particle_1_no_interaction[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
ax.scatter(x_particle_2_no_interaction[0], y_particle_2_no_interaction[0], zorder = 3, color = "#3d8c40", label = "Start particle 2")
ax.plot(x_particle_1_no_interaction, y_particle_1_no_interaction, color = "#ee5d5f", label = "Particle 1")
ax.plot(x_particle_2_no_interaction, y_particle_2_no_interaction, color = "#01377d", label = "Particle 2")
#ax.set_title("Two particle motion in x-y plane with no interactions")
ax.set_xlabel('x [$\mu$m]', fontsize = 18)
ax.set_ylabel('y [$\mu$m]', fontsize = 18)
plt.grid(color = "#e7e7e7", linestyle = "dashed")
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.axis("equal")
plt.legend(fontsize = 14, loc = "upper right")
#plt.savefig("../figures/two_particles_no_interactions_x_y.pdf")
plt.show()
##-------------------------------------------------------------------


##-------------------------------------------------------------------
#plots the two particles with interactions 
fig, ay = plt.subplots(figsize = (10,7))
ay.scatter(x_particle_1_interaction[0], y_particle_1_interaction[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
ay.scatter(x_particle_2_interaction[0], y_particle_2_interaction[0], zorder = 3, color = "#3d8c40", label = "Start particle 2")
ay.plot(x_particle_1_interaction, y_particle_1_interaction, color = "#ee5d5f", label = "Particle 1")
ay.plot(x_particle_2_interaction, y_particle_2_interaction, color = "#01377d", label = "Particle 2")
#ay.set_title("Two particle motion in x-y plane with interactions")
ay.set_xlabel('x [$\mu$m]', fontsize = 22)
ay.set_ylabel('y [$\mu$m]n', fontsize = 22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(color = "#e7e7e7", linestyle = "dashed")
plt.legend(fontsize = 15, loc = "upper right")
plt.axis("equal")
#plt.savefig("../figures/two_particles_with_interactions_x_y.pdf")
plt.show()
##-------------------------------------------------------------------



