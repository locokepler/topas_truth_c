import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numeric import full

full_data = np.loadtxt("o", delimiter=',')
relevent = full_data[:,2]
# print(relevent[0:10])

a = relevent[relevent >= 0]

# print(a[0:10])

# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

plt.hist(a, bins=80, density=True, rwidth=1., color='green')
# plt.plot(x_axis, bin_size)

plt.xlabel('Miss Distance (cm)')
plt.ylabel('Fraction of Occurances')
plt.title('Historgram of line of responce miss distance for in-patient scatters')
# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
plt.grid(False)
plt.show()



# x[x >= 0]