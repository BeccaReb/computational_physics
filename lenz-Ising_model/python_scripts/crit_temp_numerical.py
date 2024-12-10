# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats


# temperatures = [2.27, 2.28, 2.28, 2.29]
# L = [1./40., 1./60., 1./80., 1./100.]


# # slope, intercept, res, p, se = stats.linregress(L, temperatures)

# # plt.plot(L, temperatures, 'o', label='original data')
# # plt.plot(L, res.intercept + res.slope*L, 'r', label='fitted line')
# # plt.legend()
# # plt.show()


# # # linregress(temperatures, y=None)


# x = L
# y =temperatures

# m, b = np.polyfit(x, y, 1)

# plt.plot(x, y, 'yo', x, m*x+b, '--k')
# plt.show()


import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"xtick.labelsize": 13, "ytick.labelsize": 13})

x = [0.025,0.01666667,0.0125,0.01]
y = [2.28,2.29,2.29,2.30] #use when working with susceptibility 
#y = [2.28,2.28,2.29,2.27] #use when working with heat capacity

coef = np.polyfit(x,y,1)
poly1d_fn = np.poly1d(coef)
print(poly1d_fn) 

plt.scatter(x[0],y[0], color='#377eb8', label = "L = 40")
plt.scatter(x[1],y[1], color='#4daf4a', label = "L = 60")
plt.scatter(x[2],y[2], color='#e41a1c', label = "L = 80")
plt.scatter(x[3],y[3], color='#984ea3', label = "L = 100")
plt.plot(x, poly1d_fn(x), color = "dimgray")
plt.xlabel("1/L", fontsize = 14)
plt.ylabel(r"$T_c$ [ J / $k_B$]", fontsize = 14)
plt.legend(fontsize = 13)
#plt.savefig("../figures/heat_capacity_crit_temp.pdf") #use when working with heat capacity
#plt.savefig("../figures/susceptibility_crit_temp.pdf") #use when working with susceptibility 
plt.show()