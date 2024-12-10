
import numpy as np
import matplotlib.pyplot as plt

#no particle interaction case
data_particle_1_no_interaction =np.genfromtxt("../textfiles/particle_1_4000_no_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)
data_particle_2_no_interaction =np.genfromtxt("../textfiles/particle_2_4000_no_interaction_rk.txt", dtype="float", delimiter=",", skip_header=0)


time = data_particle_1_no_interaction[:,0] #time is same for all txt files!

#no particle interaction case
x_particle_1_no_interaction = data_particle_1_no_interaction[:,1]
z_particle_1_no_interaction = data_particle_1_no_interaction[:,3]
x_particle_2_no_interaction = data_particle_2_no_interaction[:,1]
z_particle_2_no_interaction = data_particle_2_no_interaction[:,3]

x_velocity_1_no_interaction = data_particle_1_no_interaction[:,4]  
z_velocity_1_no_interaction = data_particle_1_no_interaction[:,6]
x_velocity_2_no_interaction = data_particle_2_no_interaction[:,4]  
z_velocity_2_no_interaction = data_particle_2_no_interaction[:,6] 

##-------------------------------------------------------------------
#plots the two particles with no interactions (x, v_x) 
fig, ax = plt.subplots(figsize = (9,9))
ax.scatter(x_particle_1_no_interaction[0], x_velocity_1_no_interaction[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
ax.scatter(x_particle_2_no_interaction[0], x_velocity_2_no_interaction[0], zorder = 3, color = "#3d8c40", label = "Start particle 2")
ax.plot(x_particle_1_no_interaction, x_velocity_1_no_interaction, color = "#ee5d5f", label = "Particle 1")
ax.plot(x_particle_2_no_interaction, x_velocity_2_no_interaction, color = "#01377d", label = "Particle 2")
#ax.set_title("Two particles in phase space x with interactions")
ax.set_xlabel('x [$\mu$m]', fontsize = 22)
ax.set_ylabel('$v_x$ [$\mu$m / $\mu$s]', fontsize = 22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize = 16, loc = "upper right")
plt.grid(color = "#e7e7e7", linestyle = "dashed")
plt.axis("equal")
#plt.savefig("../figures/particles_no_interactions_phase_space_x.pdf")
plt.show()
##-------------------------------------------------------------------


##-------------------------------------------------------------------
#plots the two particles with no interactions (z, v_z) 
fig, ay = plt.subplots(figsize = (9,9))
ay.scatter(z_particle_1_no_interaction[0], z_velocity_1_no_interaction[0], zorder = 3, color = "#7e1662", label = "Start particle 1")
ay.scatter(z_particle_2_no_interaction[0], z_velocity_2_no_interaction[0], zorder = 3, color = "#3d8c40", label = "Start particle 2")
ay.plot(z_particle_1_no_interaction, z_velocity_1_no_interaction, color = "#ee5d5f", label = "Particle 1")
ay.plot(z_particle_2_no_interaction, z_velocity_2_no_interaction, color = "#01377d", label = "Particle 2")
#ay.set_title("Two particles in phase space x with interactions")
ay.set_xlabel('z [$\mu$m]', fontsize = 22)
ay.set_ylabel('$v_z$ [$\mu$m / $\mu$s]', fontsize = 22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize = 16, loc = "upper right")
plt.grid(color = "#e7e7e7", linestyle = "dashed")
#plt.savefig("../figures/particles_no_interactions_phase_space_z.pdf")
plt.show()
##-------------------------------------------------------------------

