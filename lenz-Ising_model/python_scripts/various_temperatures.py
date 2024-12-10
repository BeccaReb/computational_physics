
import numpy as np
import matplotlib.pyplot as plt


#load data files and arrange their data points after temperature
data_lattice_40 = np.genfromtxt("../textfiles/lattice_40.txt", delimiter = ",", dtype="float", skip_header=0)
lattice_40 = data_lattice_40[data_lattice_40[:,0].argsort()] 

data_lattice_60 = np.genfromtxt("../textfiles/lattice_60.txt", delimiter = ",", dtype="float", skip_header=0)
lattice_60 = data_lattice_60[data_lattice_60[:,0].argsort()] 

data_lattice_80 = np.genfromtxt("../textfiles/lattice_80.txt", delimiter = ",", dtype="float", skip_header=0)
lattice_80 = data_lattice_80[data_lattice_80[:,0].argsort()] 

data_lattice_100 = np.genfromtxt("../textfiles/lattice_100.txt", delimiter = ",", dtype="float", skip_header=0)
lattice_100 = data_lattice_100[data_lattice_100[:,0].argsort()] 

data_lattice_100_no_outlier = np.genfromtxt("../textfiles/lattice_100_no_outlier.txt", delimiter = ",", dtype="float", skip_header=0)
lattice_100_no_outlier = data_lattice_100_no_outlier[data_lattice_100_no_outlier[:,0].argsort()] 

#make x-axis values
x = np.arange(2.1, 2.4+0.000001, 0.01)
x_alternative = np.arange(2.11, 2.4+0.000001, 0.01) 

plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

fig, ax = plt.subplots(2,2, figsize = (12,8), constrained_layout = True, gridspec_kw = {"wspace": 0.1, "hspace": 0.12})

# epsilon
ax[0,0].plot(x, lattice_40[:,1],"--o", markersize = 4, label = "L = 40", color='#377eb8')
ax[0,0].plot(x, lattice_60[:,1],"--o", markersize = 4, label = "L = 60", color='#4daf4a')
ax[0,0].plot(x, lattice_80[:,1],"--o", markersize = 4, label = "L = 80", color='#e41a1c')
#ax[0,0].plot(x, lattice_100[:,1], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when including outlier
ax[0,0].plot(x_alternative, lattice_100_no_outlier[:,1], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when excluding outlier
ax[0,0].set_xlabel(r"T [ J / $k_B$ ]", fontsize = 14)
ax[0,0].set_ylabel(r"$\langle \epsilon \rangle$", fontsize = 14)
ax[0,0].legend(loc = "lower right")

# heat capacity
ax[0,1].plot(x, lattice_40[:,2],"--o",markersize = 4, label = "L = 40", color='#377eb8')
ax[0,1].plot(x, lattice_60[:,2],"--o",markersize = 4, label = "L = 60", color='#4daf4a')
ax[0,1].plot(x, lattice_80[:,2],"--o",markersize = 4, label = "L = 80", color='#e41a1c')
#ax[0,1].plot(x, lattice_100[:,2], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when including outlier
ax[0,1].plot(x_alternative, lattice_100_no_outlier[:,2], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when excluding outlier
ax[0,1].set_xlabel(r"T [ $J$ / $k_B$ ]", fontsize = 14)
ax[0,1].set_ylabel(r"$C_v / N $", fontsize = 14)
ax[0,1].legend()

# magnetisation
ax[1,0].plot(x, lattice_40[:,3],"--o", markersize = 4, label = "L = 40", color='#377eb8')
ax[1,0].plot(x, lattice_60[:,3],"--o", markersize = 4, label = "L = 60", color='#4daf4a')
ax[1,0].plot(x, lattice_80[:,3], "--o", markersize = 4, label = "L = 80", color='#e41a1c')
#ax[1,0].plot(x, lattice_100[:,3], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when including outlier
ax[1,0].plot(x_alternative, lattice_100_no_outlier[:,3], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when excluding outlier
ax[1,0].set_xlabel(r"T [ J / $k_B$ ]", fontsize = 14)
ax[1,0].set_ylabel(r"$\langle |m| \rangle$", fontsize = 14)
ax[1,0].legend()

# susceptibility
ax[1,1].plot(x, lattice_40[:,4], "--o", markersize = 4, label = "L = 40", color='#377eb8')
ax[1,1].plot(x, lattice_60[:,4], "--o", markersize = 4, label = "L = 60", color='#4daf4a')
ax[1,1].plot(x, lattice_80[:,4], "--o", markersize = 4, label = "L = 80", color='#e41a1c')
ax[1,1].plot(x, lattice_100[:,4], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when including outlier
ax[1,1].plot(x_alternative, lattice_100_no_outlier[:,4], "--o", markersize = 4, label = "L = 100", color='#984ea3') #use when excluding outlier
#ax[1,1].set_xlabel(r"T [ J / $k_B$ ]", fontsize = 14)
ax[1,1].set_ylabel(r"$\chi / N $", fontsize = 14)
ax[1,1].legend()


#plt.savefig("../figures/phase_transition_scatter_plot_no_outlier.pdf")
plt.show()