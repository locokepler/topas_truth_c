from typing import ChainMap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy.core.numeric import full
from numpy.lib.stride_tricks import broadcast_to
from matplotlib.ticker import MaxNLocator

full_data = np.loadtxt("30_cm_0_deg_10mil", delimiter=',')
relevent = full_data[:,1]
rel_energy = full_data[:,2]
one_scatters = full_data[full_data[:,1] == 1]
one_scat_energy = one_scatters[:,2]
rel_bright_int = full_data[:,4]
rel_bright_eng = full_data[:,3]
rel_dist1 = full_data[:,5]
rel_dist2 = full_data[:,6]
rel_dist3 = full_data[:,7]
compton_eneg = full_data[:,8]
compton_angl = full_data[:,9]
# print(relevent[0:10])

a = relevent[relevent >= 0]
b = rel_energy[rel_energy >= 0]
b = b[b != 511.]
c = one_scat_energy[one_scat_energy >= 0]
bright_int = rel_bright_int[rel_bright_int >= 0]
bright_eng = rel_bright_eng[rel_bright_eng >= 0]
rel_dist1_real = rel_dist1[rel_dist1 >= 0]
rel_dist2_real = rel_dist2[rel_dist2 >= 0]
rel_dist3_real = rel_dist3[rel_dist3 >= 0]
compton_eneg_real = compton_eneg[compton_eneg >= 0]
compton_angl_real = compton_angl[compton_angl >= 0]


# print(a[0:10])
# given_bins = np.logspace(np.log10(0.1), np.log10(100), 100)

escape_bins = 100
bright_bins = 100

# hist, bins = np.histogram(relevent, bins=given_bins, range=(0,100))
hist, bins = np.histogram(a, bins=20, range=(0,20), density=True)
escp_eng, bins2 = np.histogram(b, bins=escape_bins, range=(0, 511))
scat_1_escp_eng, bins3 = np.histogram(c, bins=escape_bins, range=(0,511))
escp_eng = (escp_eng / (len(b))) / (511. / escape_bins)
scat_1_escp_eng = (scat_1_escp_eng / (len(b))) / (511. / escape_bins)
hist2less3 = escp_eng - scat_1_escp_eng
brightest, brightbin = np.histogram(bright_int, bins=6, range=(1,7), density=True)
eng_bright, eng_bright_bin = np.histogram(bright_eng, bins=bright_bins, range=(0,340.666666))
first_scat_dist, scat_bins = np.histogram(rel_dist1_real, bins = 30, range=(0,30))
second_scat_dist, empty = np.histogram(rel_dist2_real, bins=30, range=(0,30))
third_scat_dist, empty = np.histogram(rel_dist3_real, bins=30, range=(0,30))

h_compt, y_edges, x_edges = np.histogram2d(compton_eneg_real, compton_angl_real, bins=1000)


# def compton_hist(angles, energies, bins, range):
# 	min_end = min(range)
# 	max_end = max(range)
# 	min_bin



hist = 100. * hist
escp_eng = 100. * escp_eng
scat_1_escp_eng = 100. * scat_1_escp_eng
brightest = 100. * brightest
eng_bright = 100. * ((eng_bright/ len(bright_eng)) / (340.66666 / bright_bins))
first_scat_dist = (100.*first_scat_dist) / len(rel_dist1)
second_scat_dist = (100.*second_scat_dist) / len(rel_dist1)
third_scat_dist = (100.*third_scat_dist) / len(rel_dist1)


# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

# hist, bins = np.histogram(a)
fig, ax = plt.subplots()
fig, ay = plt.subplots()
# fig, az = plt.subplots()
fig, aw = plt.subplots()
fig, an = plt.subplots()
fig, scat_dist = plt.subplots()
# fig, compton = plt.subplots()
fig, compt_2d = plt.subplots()
fig, compt_slice = plt.subplots()

# ax.scatter(hist, bins[:(len(bins) - 1)])
# ax.hist(a, bins=80, density=True, rwidth=.8, color='grey', range= (0,79), log=False)
ax.plot(bins[:-1], hist, 'Dk', label = 'Number of scatters before escape')
# ay.vlines(170.333, 0, 6.8, label='kinematic limit of one scatter', )
ay.plot(bins3[:-1], scat_1_escp_eng, '--b', mfc='none', label = 'energies of escaped gammas\nthat only scattered once')
ay.plot(bins2[:-1], escp_eng, '-k', label = 'energies of all gammas\nafter escaping')

# az.plot(bins2[:-1], hist2less3, 'ok', mfc='none')
aw.plot(brightbin[:-1], brightest, 'dk', label = 'highest energy scatter')
an.plot(eng_bright_bin[:-1], eng_bright, '-k', label = 'scatter energy')
# plt.plot(x_axis, bin_size)

scat_dist.plot(scat_bins[:-1], first_scat_dist, '-k', label='N1')
scat_dist.plot(scat_bins[:-1], second_scat_dist, '--r', label='N2')
scat_dist.plot(scat_bins[:-1], third_scat_dist, ':b', label='N3')

# compton.plot(compton_angl_real, compton_eneg_real, 'ok', markerfacecolor=(.5,.5,.5,.5))

# compt_2d.imshow(np.log(np.rot90(h_compt)), cmap='hot')
X, Y = np.meshgrid(x_edges, y_edges)
compt_2d.pcolormesh(X, Y, np.log(h_compt), cmap='hot')



compt_slice.plot(X[0,165:245], ((h_compt[165:245,:])[:,218]), label='slice at ')
print((Y[218,0]))
compt_slice.set_xlabel('angle')
compt_slice.set_ylabel('counts')
compt_slice.set_ylim(bottom=0)
#ay.set_xscale('log')
# X[180:250],
font_scale = 16

ax.set_xlabel('Scatters before escaping', fontsize=font_scale)
ax.set_ylabel('% of Occurances', fontsize=font_scale)
# ax.set_title('Number of scatters before a gamma escapes,\narriving into 30 cm of water at 45 degrees')
ay.set_xlabel('gamma energy at escape (keV)', fontsize=font_scale)
ay.set_ylabel('% of all gammas/keV', fontsize=font_scale)
# ay.set_title('Energy of escaping gammas,\narriving into 30 cm of water at 45 degrees')
aw.set_xlabel('highest energy scatter number', fontsize=font_scale)
aw.set_ylabel('% of all scatters', fontsize=font_scale)
# aw.set_title('Brightest scatters')
an.set_xlabel('energy deposited in highest energy scatter (keV)', fontsize=font_scale)
an.set_ylabel('% of all scatters with given energy/keV', fontsize=font_scale)
# an.set_title('Energy deposited by brightest scatter')

scat_dist.set_xlabel('cm into water', fontsize=font_scale)
scat_dist.set_ylabel('% of events/cm', fontsize=font_scale)

compt_2d.set_xlabel('scatter angle', fontsize=font_scale)
compt_2d.set_ylabel('deposited energy', fontsize=font_scale)

ax.set_ylim(bottom=0.)
ay.set_ylim(bottom=0.)
# ay.set_ylim(top=6.8)
aw.set_ylim(bottom=0.)
an.set_ylim(bottom=0.)
scat_dist.set_ylim(bottom=0., top=10.)
scat_dist.set_xlim(left=0., right=29.)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True);
ax.tick_params(bottom=True, top=True, left=True, right=True)

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
# ax.legend()
ay.legend()
# aw.legend()
# an.legend()
scat_dist.legend()
plt.tick_params(bottom=True, top=True, left=True, right=True)
# plt.tick_params(which = 'minor',bottom=False, top=False, left=False, right=False)
plt.show()



# x[x >= 0]
