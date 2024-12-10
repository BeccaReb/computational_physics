import numpy as np
import matplotlib.pyplot as plt


# names of data files we will use
data_files = ("quantities_10_0.txt","quantities_10_1.txt","quantities_10_2.txt", "quantities_10_3.txt", "quantities_10_4.txt", "quantities_10_5.txt", "quantities_10_6.txt")

#load and print quantities from data files
for i in range(0,7,1):
    data = np.genfromtxt(f"../textfiles/{data_files[i]}", dtype="float", delimiter = ",", skip_header=0)
    print("\n")
    print(f"For 10^{i} cycles, the average energy per spin is ", np.round(np.average(data[:,0]), 5))
    print(f"For 10^{i} cycles, the average magnetisation per spin is ", np.round(np.average(data[:,2]), 5))
    print(f"For 10^{i} cycles, the average heat capacity is ", np.round(np.average(data[:,1]), 5))
    print(f"For 10^{i} cycles, the average susceptibility is ", np.round(np.average(data[:,3]), 5))

print("\n")

    