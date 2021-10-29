import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy.core.numeric import full
from numpy.lib.stride_tricks import broadcast_to

full_data = np.loadtxt("30_cm_45_deg", delimiter=',')
relevent = full_data[:,1]
rel_energy = full_data[:,2]
one_scatters = full_data[full_data[:,1] == 1]
one_scat_energy = one_scatters[:,2]
rel_bright_int = full_data[:,4]
rel_bright_eng = full_data[:,3]
# print(relevent[0:10])

a = relevent[relevent >= 0]
b = rel_energy[rel_energy >= 0]
c = one_scat_energy[one_scat_energy >= 0]
bright_int = rel_bright_int[rel_bright_int >= 0]
bright_eng = rel_bright_eng[rel_bright_eng >= 0]

# print(a[0:10])
# given_bins = np.logspace(np.log10(0.1), np.log10(100), 100)

# hist, bins = np.histogram(relevent, bins=given_bins, range=(0,100))
hist, bins = np.histogram(a, bins=20, range=(0,20), normed=True)
hist2, bins2 = np.histogram(b, bins=40, range=(0, 511), normed=False)
hist3, bins3 = np.histogram(c, bins=40, range=(0,511), normed=False)
hist2 = hist2 / len(b)
hist3 = hist3 / len(b)
hist2less3 = hist2 - hist3
brightest, brightbin = np.histogram(bright_int, bins=6, range=(1,7), normed=True)
eng_bright, eng_bright_bin = np.histogram(bright_eng, bins=40, range=(0,340.666666), normed=True)

hist = 100 * hist
hist2 = 100 * hist2
hist3 = 100 * hist3
brightest = 100 * brightest
eng_bright = 100 * eng_bright


# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

# hist, bins = np.histogram(a)
fig, ax = plt.subplots()
fig, ay = plt.subplots()
fig, az = plt.subplots()
fig, aw = plt.subplots()
fig, an = plt.subplots()

# ax.scatter(hist, bins[:(len(bins) - 1)])
# ax.hist(a, bins=80, density=True, rwidth=.8, color='grey', range= (0,79), log=False)
ax.plot(bins[:-1], hist, 'Dk', label = 'Number of scatters before escape of 511 kev\ngamma entering 30 cm of water at 45 degrees')
ay.vlines(170.333, 0, 6.8, label='kinematic limit of one scatter', )
ay.plot(bins3[:-1], hist3, 'ok', mfc='none', label = 'energies for gammas at 511 keV\nthat escape after one scatter')
ay.plot(bins2[:-1], hist2, 'Dk-', label = 'energies for all gammas after\nescaping after entering a 30 cm\nthick block of water')

# az.plot(bins2[:-1], hist2less3, 'ok', mfc='none')
aw.plot(brightbin[:-1], brightest, 'dk', label = 'brightest scatter by a\ngamma at 511 keV entering\nwater normal to the surface')
an.plot(eng_bright_bin[:-1], eng_bright, 'dk', label = 'highest energy scatter of a\ngamma at 511 keV entering\nwater normal to the surface')
# plt.plot(x_axis, bin_size)
#ay.set_xscale('log')

ax.set_xlabel('Scatters before escaping')
ax.set_ylabel('Percent of Occurances')
# ax.set_title('Number of scatters before a gamma escapes,\narriving into 30 cm of water at 45 degrees')
ay.set_xlabel('gamma energy at escape (keV)')
ay.set_ylabel('Percent of all gammas')
# ay.set_title('Energy of escaping gammas,\narriving into 30 cm of water at 45 degrees')
aw.set_xlabel('brightest scatter number')
aw.set_ylabel('Percent of all scatters')
# aw.set_title('Brightest scatters')
an.set_xlabel('energy deposited in brightest scatter (keV)')
an.set_ylabel('Percent of all scatters with given energy')
# an.set_title('Energy deposited by brightest scatter')

ax.set_ylim(bottom=0.)
ay.set_ylim(bottom=0.)
ay.set_ylim(top=6.8)
aw.set_ylim(bottom=0.)
an.set_ylim(bottom=0.)


# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
plt.grid(False)
#from matplotlib.ticker import AutoMinorLocator
#ax.xaxis.set_minor_locator(AutoMinorLocator(10))
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ay.yaxis.set_minor_locator(AutoMinorLocator(10))
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
#ay.xaxis.set_minor_locator(locmin)
#ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.legend()
ay.legend()
aw.legend()
an.legend()
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(which = 'minor',bottom=False, top=False, left=False, right=False)
plt.show()



# x[x >= 0]
