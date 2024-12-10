
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl

#load data files
#this is for just cycles
eps_T10_u = np.genfromtxt("../textfiles/cycles_T10_u.txt", dtype="float", skip_header=0)
eps_T10_o = np.genfromtxt("../textfiles/cycles_T10_o.txt", dtype="float", skip_header=0) 
eps_T24_u = np.genfromtxt("../textfiles/cycles_T24_u.txt", dtype="float", skip_header=0) 
eps_T24_o = np.genfromtxt("../textfiles/cycles_T24_o.txt", dtype="float", skip_header=0) 

#this if for both samples and cycles
mean_eps_T10_u_data = np.genfromtxt("../textfiles/cycles_samples_T10_u.txt", dtype="float", skip_header=0)
mean_eps_T10_o_data = np.genfromtxt("../textfiles/cycles_samples_T10_o.txt", dtype="float", skip_header=0)
mean_eps_T24_u_data = np.genfromtxt("../textfiles/cycles_samples_T24_u.txt", dtype="float", skip_header=0)
mean_eps_T24_o_data = np.genfromtxt("../textfiles/cycles_samples_T24_o.txt", dtype="float", skip_header=0)


n_MC_cycles = np.linspace(1, len(mean_eps_T10_u_data), len(mean_eps_T10_u_data)) #we have run all variations with 10^5 cycles, so their lengths should be the same and it doesnt matter which one we choose in this line

#calculate averages
mean_eps_T10_u = []
energy_sum_T10_u = 0
i = 0    
while i < len(n_MC_cycles):
    energy_sum_T10_u = (np.sum(mean_eps_T10_u_data[0:i+1]))
    mean_eps_T10_u.append(energy_sum_T10_u/n_MC_cycles[i])
    i+=1
    
mean_eps_T10_o = []
energy_sum_T10_o = 0
i = 0    
while i < len(n_MC_cycles):
    energy_sum_T10_o = (np.sum(mean_eps_T10_o_data[0:i+1]))
    mean_eps_T10_o.append(energy_sum_T10_o/n_MC_cycles[i])
    i+=1

mean_eps_T24_u = []
energy_sum_T24_u = 0
i = 0    
while i < len(n_MC_cycles):
    energy_sum_T24_u = (np.sum(mean_eps_T24_u_data[0:i+1]))
    mean_eps_T24_u.append(energy_sum_T24_u/n_MC_cycles[i])
    i+=1

mean_eps_T24_o = []
energy_sum_T24_o = 0
i = 0    
while i < len(n_MC_cycles):
    energy_sum_T24_o = (np.sum(mean_eps_T24_o_data[0:i+1]))
    mean_eps_T24_o.append(energy_sum_T24_o/n_MC_cycles[i])
    i+=1

# plot figure   
plt.figure(figsize = (24, 16))

plt.plot(np.log10(n_MC_cycles), eps_T10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0, label = r"$\epsilon$, T=1.0 $J/k_{B}$, unordered", zorder = 8)
plt.plot(np.log10(n_MC_cycles), eps_T10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0, label = r"$\epsilon$, T=1.0 $J/k_{B}$, ordered")
plt.plot(np.log10(n_MC_cycles), eps_T24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0, label = r"$\epsilon$, T=2.4 $J/k_{B}$, unordered", zorder = 7)
plt.plot(np.log10(n_MC_cycles), eps_T24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0, label = r"$\epsilon$, T=2.4 $J/k_{B}$, ordered")

plt.plot(np.log10(n_MC_cycles), mean_eps_T10_u, '-', linewidth=2.0, color='#377eb8', label=r'$\langle \epsilon \rangle$, $T=1.0$ $J/k_{B}$, unordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T10_o, '-', linewidth=2.0, color='#4daf4a', label=r'$\langle \epsilon \rangle$, $T=1.0$ $J/k_{B}$, ordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c', label=r'$\langle \epsilon \rangle$, $T=2.4$ $J/k_{B}$, unordered')
plt.plot(np.log10(n_MC_cycles), mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3', label=r'$\langle \epsilon \rangle$, $T=2.4$ $J/k_{B}$, ordered')

plt.xticks(size = 26)
plt.yticks(size = 26)
plt.xlabel("Cycles",fontsize = 30)
plt.ylabel(r"[$\;$J$\;$]",fontsize = 30)
plt.legend(fontsize = 22, loc = "upper right")

#plt.savefig("../figures/burn_in_time.pdf")
plt.show()

################## Probability function ##################
mpl.rcParams['font.size'] = 13
fig, ax = plt.subplots(2, figsize = (6, 10))
hist_10 = ax[0].hist(mean_eps_T10_u_data[10**2:10**5], 100, label = "T=1.0 $J/k_B$")
hist_24 = ax[1].hist(mean_eps_T24_u_data[10**2:10**5], 100, label = "T=2.4 $J/k_B$")

def gaussian(x, a, b, c, d):
    """
    a: the amplitude of the gaussian.
    b: the mean of the gaussian. The position of the peak on the x-axis.
    c: the standard deviation.
    d: the constant term, y-value of the baseline.
    """
    return a*np.exp(-((x - b)**2) / (2*c**2)) + d

x = np.linspace(-1.8, -0.8, 100)
popt, pcov = curve_fit(gaussian, x, hist_24[0], p0 = [2500, -1.2, 0.1, 0])
ax[1].plot(x, gaussian(x, *popt))

fig.supylabel("Counts")
fig.supxlabel(r"$\langle \epsilon \rangle$")
ax[0].legend(); ax[1].legend()
#plt.savefig("../probability_function.pdf")
plt.show()



   