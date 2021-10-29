import numpy as np
import matplotlib.pyplot as plt
import matplotlib

full_data = np.loadtxt("1%_cut", delimiter=',')
no_in_pat = full_data#[full_data[:,1] == 0]
relevent = no_in_pat[:,2]
# print(relevent[0:10])

a = relevent[relevent >= 0]

# print(a[0:10])
# given_bins = np.logspace(np.log10(0.1), np.log10(100), 100)

# hist, bins = np.histogram(relevent, bins=given_bins, range=(0,100))
hist, bins = np.histogram(relevent, bins=60, range=(0,60))

hist2, bins2 = np.histogram(relevent, bins= np.logspace(np.log10(0.001),np.log10(100.), 40), range=(0, 100))

# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

# hist, bins = np.histogram(a)
fig, ax = plt.subplots()
fig, ay = plt.subplots()

# ax.scatter(hist, bins[:(len(bins) - 1)])
# ax.hist(a, bins=80, density=True, rwidth=.8, color='grey', range= (0,79), log=False)
ax.plot(bins[:-1], hist, 'Dk')
ay.plot(bins2[:-1], hist2, 'Dk')
# plt.plot(x_axis, bin_size)
ay.set_xscale('log')

ax.set_xlabel('Miss Distance (cm)')
ax.set_ylabel('Number of Occurances')
ax.set_title('Histogram of line of responce miss distance for reconstructed scatters\n1% energy cut')
ay.set_xlabel('Miss Distance (cm)')
ay.set_ylabel('Number of Occurances')
ay.set_title('Histogram of line of responce miss distance for reconstructed scatters')

# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
plt.grid(False)
from matplotlib.ticker import AutoMinorLocator
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ay.yaxis.set_minor_locator(AutoMinorLocator(10))
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
ay.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(which = 'minor',bottom=True, top=True, left=True, right=True)
plt.show()



# x[x >= 0]
