import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numeric import full

full_data = np.loadtxt("20_cm_head_on", delimiter=',')
#no_in_pat = full_data[full_data[:,1] == 0]
relevent = full_data[:,1] #no_in_pat[:,2]
# print(relevent[0:10])

a = relevent[relevent >= 0]

# print(a[0:10])

# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

# hist, bins = np.histogram(a)
fig, ax = plt.subplots()
# ax.scatter(hist, bins[:(len(bins) - 1)])
ax.hist(a, bins=80, density=True, rwidth=.8, color='grey', range= (0,79), log=True)
# plt.plot(x_axis, bin_size)

plt.xlabel('Miss Distance (cm)')
plt.ylabel('Fraction of Occurances')
plt.title('Historgram of line of responce miss distance for scatters')
# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
plt.grid(False)
from matplotlib.ticker import AutoMinorLocator
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(which = 'minor',bottom=True, top=True, left=True, right=True)
plt.show()



# x[x >= 0]
