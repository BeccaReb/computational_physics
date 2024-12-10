
import numpy as np
import matplotlib.pyplot as plt


#this is for just cycles
eps_T10_u = np.genfromtxt("../textfiles/cycles_T10_u.txt", dtype="float", skip_header=0)
eps_T10_o = np.genfromtxt("../textfiles/cycles_T10_o.txt", dtype="float", skip_header=0) 
eps_T24_u = np.genfromtxt("../textfiles/cycles_T24_u.txt", dtype="float", skip_header=0) 
eps_T24_o = np.genfromtxt("../textfiles/cycles_T24_o.txt", dtype="float", skip_header=0) 

#this if for both samples and cycles
mean_eps_T10_u_data = np.genfromtxt("../textfiles/cycles_samples_T10_u_test.txt", dtype="float", skip_header=0)
mean_eps_T10_o_data = np.genfromtxt("../textfiles/cycles_samples_T10_o_test.txt", dtype="float", skip_header=0)
mean_eps_T24_u_data = np.genfromtxt("../textfiles/cycles_samples_T24_u_test.txt", dtype="float", skip_header=0)
mean_eps_T24_o_data = np.genfromtxt("../textfiles/cycles_samples_T24_o_test.txt", dtype="float", skip_header=0)


n_MC_cycles = np.linspace(1, len(mean_eps_T10_u_data), len(mean_eps_T10_u_data)) #we have run all variations with 10^5 cycles, so their lengths should be the same and it doesnt matter which one we choose in this line

print(mean_eps_T10_u_data[0:20])

# mean_eps_T10_u = []
# energy_sum_T10_u = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T10_u = (np.sum(mean_eps_T10_u_data[0:i+1]))
#     mean_eps_T10_u.append(energy_sum_T10_u/n_MC_cycles[i])
#     i+=1

    
# mean_eps_T10_o = []
# energy_sum_T10_o = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T10_o = (np.sum(mean_eps_T10_o_data[0:i+1]))
#     mean_eps_T10_o.append(energy_sum_T10_o/n_MC_cycles[i])
#     i+=1

# mean_eps_T24_u = []
# energy_sum_T24_u = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T24_u = (np.sum(mean_eps_T24_u_data[0:i+1]))
#     mean_eps_T24_u.append(energy_sum_T24_u/n_MC_cycles[i])
#     i+=1

# mean_eps_T24_o = []
# energy_sum_T24_o = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T24_o = (np.sum(mean_eps_T24_o_data[0:i+1]))
#     mean_eps_T24_o.append(energy_sum_T24_o/n_MC_cycles[i])
#     i+=1

# print(mean_eps_T10_u[0:20])
    
# mean_eps_T10_o = []
# energy_sum_T10_o = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T10_o = (mean_eps_T10_o_data[i]/(n_MC_cycles[i]))
#     mean_eps_T10_o.append(energy_sum_T10_o)
#     i+=1
    
# mean_eps_T24_u = []
# energy_sum_T24_u = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T24_u = (mean_eps_T24_u_data[i]/(n_MC_cycles[i]))
#     mean_eps_T24_u.append(energy_sum_T24_u)
#     i+=1

# mean_eps_T24_o = []
# energy_sum_T24_o = 0
# i = 0    
# while i < len(n_MC_cycles):
#     energy_sum_T24_o = (mean_eps_T24_o_data[i]/(n_MC_cycles[i]))
#     mean_eps_T24_o.append(energy_sum_T24_o)
#     i+=1

# print((n_MC_cycles))
    
plt.figure()

plt.plot(np.log10(n_MC_cycles), eps_T10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
plt.plot(np.log10(n_MC_cycles), eps_T10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0)
plt.plot(np.log10(n_MC_cycles), eps_T24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
plt.plot(np.log10(n_MC_cycles), eps_T24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0)

plt.plot(np.log10(n_MC_cycles), mean_eps_T10_u_data, '-', linewidth=2.0, color='#377eb8', label='$T=1.0$ $J/k_{B}$, unordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T10_o_data, '-', linewidth=2.0, color='#4daf4a', label='$T=1.0$ $J/k_{B}$, ordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T24_u_data, '-', linewidth=2.0, color='#e41a1c', label='$T=2.4$ $J/k_{B}$, unordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T24_o_data, '-', linewidth=2.0, color='#984ea3', label='$T=2.4$ $J/k_{B}$, ordered')
plt.legend()

# plt.savefig("../figures/burn.pdf")
plt.show()

#out by a factor of 4???


#------------------------------------------------------------------------
#maybe use later:
    
# print(len(energy_vec))
#print(energy_vec)


# # print(energy_vec)
# print(len(energy_vec))
# print(energy_vec[-1])


# E_vec = []
# E = 0
# j = 0
# while j < len(n_MC_cycles):
#     E = (energy_vec[j]) / n_MC_cycles[j] 
#     E_vec.append(E) 
#     j+=1 

# print(len(E_vec))


# plt.hist(eps_T10_u, 1000)
# plt.show()


   