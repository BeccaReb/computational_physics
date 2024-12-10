import numpy as np
import matplotlib.pyplot as plt


#loading analytical datafiles
analytical_4000 =np.genfromtxt("../textfiles/single_particle_1_4000_analytical.txt", dtype="float", delimiter=",", skip_header=0)
analytical_8000 =np.genfromtxt("../textfiles/single_particle_1_8000_analytical.txt", dtype="float", delimiter=",", skip_header=0)
analytical_16000 =np.genfromtxt("../textfiles/single_particle_1_16000_analytical.txt", dtype="float", delimiter=",", skip_header=0)
analytical_32000 =np.genfromtxt("../textfiles/single_particle_1_32000_analytical.txt", dtype="float", delimiter=",", skip_header=0)

#loading forward euler data files
euler_4000 =np.genfromtxt("../textfiles/single_particle_1_4000_euler.txt", dtype="float", delimiter=",", skip_header=0)
euler_8000 =np.genfromtxt("../textfiles/single_particle_1_8000_euler.txt", dtype="float", delimiter=",", skip_header=0)
euler_16000 =np.genfromtxt("../textfiles/single_particle_1_16000_euler.txt", dtype="float", delimiter=",", skip_header=0)
euler_32000 =np.genfromtxt("../textfiles/single_particle_1_32000_euler.txt", dtype="float", delimiter=",", skip_header=0)

#loading RK4 data files
rk_4000 =np.genfromtxt("../textfiles/single_particle_1_4000_rk.txt", dtype="float", delimiter=",", skip_header=0)
rk_8000 =np.genfromtxt("../textfiles/single_particle_1_8000_rk.txt", dtype="float", delimiter=",", skip_header=0)
rk_16000 =np.genfromtxt("../textfiles/single_particle_1_16000_rk.txt", dtype="float", delimiter=",", skip_header=0)
rk_32000 =np.genfromtxt("../textfiles/single_particle_1_32000_rk.txt", dtype="float", delimiter=",", skip_header=0)

#calculating r for analytical data
r_analytical_4000 = np.sqrt((analytical_4000[:,1])**2 + (analytical_4000[:,2])**2 + (analytical_4000[:,3])**2)
r_analytical_8000 = np.sqrt((analytical_8000[:,1])**2 + (analytical_8000[:,2])**2 + (analytical_8000[:,3])**2)
r_analytical_16000 = np.sqrt((analytical_16000[:,1])**2 + (analytical_16000[:,2])**2 + (analytical_16000[:,3])**2)
r_analytical_32000 = np.sqrt((analytical_32000[:,1])**2 + (analytical_32000[:,2])**2 + (analytical_32000[:,3])**2)

#calculating r for euler data
r_euler_4000 = np.sqrt((euler_4000[:,1])**2 + (euler_4000[:,2])**2 + (euler_4000[:,3])**2)
r_euler_8000 = np.sqrt((euler_8000[:,1])**2 + (euler_8000[:,2])**2 + (euler_8000[:,3])**2)
r_euler_16000 = np.sqrt((euler_16000[:,1])**2 + (euler_16000[:,2])**2 + (euler_16000[:,3])**2)
r_euler_32000 = np.sqrt((euler_32000[:,1])**2 + (euler_32000[:,2])**2 + (euler_32000[:,3])**2)

#calculating r for RK data
r_rk_4000 = np.sqrt((rk_4000[:,1])**2 + (rk_4000[:,2])**2 + (rk_4000[:,3])**2)
r_rk_8000 = np.sqrt((rk_8000[:,1])**2 + (rk_8000[:,2])**2 + (rk_8000[:,3])**2)
r_rk_16000 = np.sqrt((rk_16000[:,1])**2 + (rk_16000[:,2])**2 + (rk_16000[:,3])**2)
r_rk_32000 = np.sqrt((rk_32000[:,1])**2 + (rk_32000[:,2])**2 + (rk_32000[:,3])**2)


#time arrays for different time steps
time_4000 = analytical_4000[:,0]
time_8000 = analytical_8000[:,0]
time_16000 = analytical_16000[:,0]
time_32000 = analytical_32000[:,0]

#calculating relative error RK4
error_rk_4000 = np.abs(r_analytical_4000-r_rk_4000)/r_analytical_4000
error_rk_8000 = np.abs(r_analytical_8000-r_rk_8000)/r_analytical_8000
error_rk_16000 = np.abs(r_analytical_16000-r_rk_16000)/r_analytical_16000
error_rk_32000 = np.abs(r_analytical_32000-r_rk_32000)/r_analytical_32000

#calculating relative error euler
error_euler_4000 = np.abs(r_analytical_4000-r_euler_4000)/r_analytical_4000
error_euler_8000 = np.abs(r_analytical_8000-r_euler_8000)/r_analytical_8000
error_euler_16000 = np.abs(r_analytical_16000-r_euler_16000)/r_analytical_16000
error_euler_32000 = np.abs(r_analytical_32000-r_euler_32000)/r_analytical_32000


#--------------------------------------------------------------------------------------
#Plotting relative error 
 
#plotting relative error RK4
fig, ax = plt.subplots(figsize = (9,7))
ax.plot(time_4000, error_rk_4000, color = "#fb9062", label = "4000 steps")
ax.plot(time_8000, error_rk_8000, color = "#ee5d5f", label = "8000 steps")
ax.plot(time_16000, error_rk_16000, color = "#ce4993", label = "16000 steps")
ax.plot(time_32000, error_rk_32000, color = "#6a0d83", label = "32000 steps")
ax.set_xlabel("time [$\mu$s]", fontsize = 18)
ax.set_ylabel("relative error", fontsize = 18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#ax.set_title(r"Relative error in $\mathbf{r_i}$ for Runge-Kutta 4 method")
plt.grid(color = "#e7e7e7", linestyle = "dashed")
plt.legend(fontsize = 16)
#plt.savefig("../figures/relative_error_rk4.pdf")
plt.show()


#plotting relative error euler
fig, ay = plt.subplots(figsize = (9,7))
ay.plot(time_4000, error_euler_4000, color = "#01377d", label = "4000 steps")
ay.plot(time_8000, error_euler_8000, color = "#009dd1", label = "8000 steps")
ay.plot(time_16000, error_euler_16000, color = "#37782c", label = "16000 steps")
ay.plot(time_32000, error_euler_32000, color = "#9fd283", label = "32000 steps")
ay.set_xlabel("time [$\mu$s]", fontsize = 18)
ay.set_ylabel("relative error", fontsize = 18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#ay.set_title(r"Relative error in $\mathbf{r_i}$ for Forward Euler method")
plt.legend(fontsize = 16)
plt.grid(color = "#e7e7e7", linestyle = "dashed")
#plt.savefig("../figures/relative_error_euler.pdf")
plt.show()

#--------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#Error convergence rate Runge-Kutta:
delta_max_2_rk = np.max(np.abs(r_analytical_8000-r_rk_8000))
delta_max_3_rk = np.max(np.abs(r_analytical_16000-r_rk_16000))
delta_max_4_rk = np.max(np.abs(r_analytical_32000-r_rk_32000))

delta_max_2_minus_1_rk = np.max(np.abs(r_analytical_4000-r_rk_4000)) 
delta_max_3_minus_1_rk = delta_max_2_rk
delta_max_4_minus_1_rk = delta_max_3_rk

h_2_rk = 50/8000
h_3_rk = 50/16000
h_4_rk = 50/32000

h_2_minus_1_rk = (h_2_rk)*2
h_3_minus_1_rk = h_2_rk
h_4_minus_1_rk = h_3_rk

r_err_2_rk = (np.log10(delta_max_2_rk / delta_max_2_minus_1_rk))/(np.log10(h_2_rk / h_2_minus_1_rk))
r_err_3_rk = (np.log10(delta_max_3_rk / delta_max_3_minus_1_rk))/(np.log10(h_3_rk / h_3_minus_1_rk))
r_err_4_rk = (np.log10(delta_max_4_rk / delta_max_4_minus_1_rk))/(np.log10(h_4_rk / h_4_minus_1_rk))

total_error_convergence_rate_rk = (1/3)*(r_err_2_rk + r_err_3_rk + r_err_4_rk)

print("The error convergence rate for the Runge-Kutta method is",total_error_convergence_rate_rk)

#---------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------
#Error convergence rate Euler:
delta_max_2_euler = np.max(np.abs(r_analytical_8000-r_euler_8000))
delta_max_3_euler = np.max(np.abs(r_analytical_16000-r_euler_16000))
delta_max_4_euler = np.max(np.abs(r_analytical_32000-r_euler_32000))

delta_max_2_minus_1_euler = np.max(np.abs(r_analytical_4000-r_euler_4000)) 
delta_max_3_minus_1_euler = delta_max_2_euler
delta_max_4_minus_1_euler = delta_max_3_euler

h_2_euler = 50/8000
h_3_euler = 50/16000
h_4_euler = 50/32000

h_2_minus_1_euler = (h_2_euler)*2
h_3_minus_1_euler = h_2_euler
h_4_minus_1_euler = h_3_euler

r_err_2_euler = (np.log10(delta_max_2_euler / delta_max_2_minus_1_euler))/(np.log10(h_2_euler / h_2_minus_1_euler))
r_err_3_euler = (np.log10(delta_max_3_euler / delta_max_3_minus_1_euler))/(np.log10(h_3_euler / h_3_minus_1_euler))
r_err_4_euler = (np.log10(delta_max_4_euler / delta_max_4_minus_1_euler))/(np.log10(h_4_euler / h_4_minus_1_euler))

total_error_convergence_rate_euler = (1/3)*(r_err_2_euler + r_err_3_euler + r_err_4_euler)

print("The error convergance rate for the Forward Euler method is",total_error_convergence_rate_euler)

#--------------------------------------------------------------------------------------------------------------------- 

