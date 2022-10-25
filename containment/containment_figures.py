import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker
from matplotlib.ticker import (MaxNLocator, MultipleLocator, AutoMinorLocator)

def plot_prettier(dpi=200, fontsize=17):
	plt.rcParams['figure.dpi'] = dpi
	plt.rc("savefig", dpi=dpi)
	plt.rc('font', size=fontsize)
	plt.rc('xtick', direction='in')
	plt.rc('ytick', direction='in')
	plt.rc('xtick.major', pad=5)
	plt.rc('xtick.minor', pad=5)
	plt.rc('ytick.major', pad=5)
	plt.rc('ytick.minor', pad=5)
	plt.rc('lines', dotted_pattern = [2., 2.])
	# plt.rc('text', usetex=True)

plot_prettier()

def add_minor_ticks(plot, ticks, bot=True, tp=True, lft=True, rght=True):
	plot.tick_params(which = 'minor', bottom=bot, top=tp, left=lft, right=rght)
	plot.tick_params(bottom=True, top=True, left=True, right=True)
	if (bot or tp):
		plot.xaxis.set_minor_locator((AutoMinorLocator(int(ticks))))
	if (lft or rght):
		plot.yaxis.set_minor_locator((AutoMinorLocator(int(ticks))))

# full_data = np.loadtxt("containment/30_cm_0_deg_10mil", delimiter=',')
# full_data = np.loadtxt("containment/30_cm_0_deg", delimiter=',')
# full_data = np.loadtxt("containment/thick_0deg.tuple", delimiter=',')
# full_data = np.loadtxt("containment/30_cm_0_deg_10mil_LAB.tuple", delimiter=',')
full_data = np.loadtxt("containment/very_thick_head_on_10mil.tuple", delimiter=',')
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
hist, bins = np.histogram(a, bins=20, range=(0,20))
escp_eng, bins2 = np.histogram(b, bins=escape_bins, range=(0, 511))
escp_eng_err = np.sqrt(escp_eng)
scat_1_escp_eng, bins3 = np.histogram(c, bins=escape_bins, range=(0,511))
scat_1_escp_eng_err = np.sqrt(scat_1_escp_eng)
escp_eng = (escp_eng / (len(b))) / (511. / escape_bins)
scat_1_escp_eng = (scat_1_escp_eng / (len(b))) / (511. / escape_bins)
escp_eng_err = (escp_eng_err / (len(b))) / (511. / escape_bins)
scat_1_escp_eng_err = (scat_1_escp_eng_err / (len(b))) / (511. / escape_bins)
hist2less3 = escp_eng - scat_1_escp_eng
brightest, brightbin = np.histogram(bright_int, bins=6, range=(1,7))
brightest_err = np.sqrt(brightest)
brightest = brightest / len(bright_int)
brightest_err = brightest_err / len(bright_int)
eng_bright, eng_bright_bin = np.histogram(bright_eng, bins=bright_bins, range=(0,340.666666))
eng_bright_err = np.sqrt(eng_bright)
first_scat_dist, scat_bins = np.histogram(rel_dist1_real, bins = 30, range=(0,30))
second_scat_dist, empty = np.histogram(rel_dist2_real, bins=30, range=(0,30))
third_scat_dist, empty = np.histogram(rel_dist3_real, bins=30, range=(0,30))

h_compt, y_edges, x_edges = np.histogram2d(compton_eneg_real, compton_angl_real, bins=1000)


# def compton_hist(angles, energies, bins, range):
# 	min_end = min(range)
# 	max_end = max(range)
# 	min_bin



hist_norm = 100. * hist / len(relevent)
escp_eng = 100. * escp_eng
scat_1_escp_eng = 100. * scat_1_escp_eng
escp_eng_err = 100. * escp_eng_err
scat_1_escp_eng_err = 100. * scat_1_escp_eng_err
brightest = 100. * brightest
brightest_err = 100. * brightest_err
eng_bright = 100. * ((eng_bright/ len(bright_eng)) / (340.66666 / bright_bins))
eng_bright_err = 100. * ((eng_bright_err/ len(bright_eng)) / (340.66666 / bright_bins))
first_scat_dist = (100.*first_scat_dist) / len(rel_dist1)
second_scat_dist = (100.*second_scat_dist) / len(rel_dist1)
third_scat_dist = (100.*third_scat_dist) / len(rel_dist1)


# bin_vals = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]

# bin_size = np.histogram(a, bins=20)
# print(bin_size)

# x_axis = np.linspace(1, len(bin_size),len(bin_size))

# hist, bins = np.histogram(a)
fig_ax, ax = plt.subplots()
plt.subplots_adjust(bottom=.15)
fig_ay, ay = plt.subplots()
plt.subplots_adjust(bottom=.15)
# fig, az = plt.subplots()
fig_aw, aw = plt.subplots()
plt.subplots_adjust(bottom=.15)
fig_an, an = plt.subplots()
plt.subplots_adjust(bottom=.15)
fig_scat_dist, scat_dist = plt.subplots()
plt.subplots_adjust(bottom=.15)
# fig, compton = plt.subplots()
fig_compt_2d, compt_2d = plt.subplots()
plt.subplots_adjust(bottom=.15)
plt.subplots_adjust(left=.14)
fig_compt_slice, compt_slice = plt.subplots()
plt.subplots_adjust(bottom=.15)
plt.subplots_adjust(right=.98)
plt.subplots_adjust(left=.25)
fig_compt_slice.set_size_inches(4.0,4.8)

# ax.scatter(hist, bins[:(len(bins) - 1)])
# ax.hist(a, bins=80, density=True, rwidth=.8, color='grey', range= (0,79), log=False)
# ax.plot(bins[:-1], hist, 'Dk', label = 'Number of scatters before escape')
ax.errorbar(bins[:-1], hist_norm, yerr=np.sqrt(hist) * 100 / len(a), marker='D', linestyle="None", label = 'Number of scatters\nbefore escape')
# ay.vlines(170.333, 0, 6.8, label='kinematic limit of one scatter', )
# ay.plot(bins3[:-1], scat_1_escp_eng, '--b', mfc='none', label = 'once scattered')
ay.errorbar(bins3[:-1], scat_1_escp_eng, yerr=scat_1_escp_eng_err, marker='D', mfc='none', label = 'once scattered')
# ay.plot(bins2[:-1], escp_eng, '-k', label = 'all gammas')
ay.errorbar(bins2[:-1], escp_eng, yerr=escp_eng_err, marker='s', label = 'all gammas')

# az.plot(bins2[:-1], hist2less3, 'ok', mfc='none')
# aw.plot(brightbin[:-1], brightest, 'dk', label = 'highest energy scatter')
aw.errorbar(brightbin[:-1], brightest, yerr=brightest_err, marker='D', linestyle="None", label = 'highest energy scatter')
print('brightest scatters:')
print(brightest)

# an.plot(eng_bright_bin[:-1], eng_bright, '-k', label = 'scatter energy')
an.errorbar(eng_bright_bin[:-1], eng_bright, yerr=eng_bright_err, marker='D', label = 'scatter energy')

# plt.plot(x_axis, bin_size)

add_minor_ticks(ax, 5, bot=True,tp=True)

# ax.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# ax.tick_params(bottom=True, top=True, left=True, right=True)
# ax.xaxis.set_minor_locator((AutoMinorLocator(5)))
# ax.yaxis.set_minor_locator((AutoMinorLocator(5)))


# ay.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# ay.tick_params(bottom=True, top=True, left=True, right=True)
# ay.xaxis.set_minor_locator((AutoMinorLocator(5)))
# ay.yaxis.set_minor_locator((AutoMinorLocator(5)))

add_minor_ticks(ay, 5)

# aw.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# aw.tick_params(bottom=True, top=True, left=True, right=True)
# aw.xaxis.set_minor_locator((AutoMinorLocator(5)))
# aw.yaxis.set_minor_locator((AutoMinorLocator(5)))

add_minor_ticks(aw, 5, bot=False,tp=False)

# an.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# an.tick_params(bottom=True, top=True, left=True, right=True)
# an.xaxis.set_minor_locator((AutoMinorLocator(5)))
# an.yaxis.set_minor_locator((AutoMinorLocator(5)))

add_minor_ticks(an, 5)

# scat_dist.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# scat_dist.tick_params(bottom=True, top=True, left=True, right=True)
# scat_dist.xaxis.set_minor_locator((AutoMinorLocator(5)))
# scat_dist.yaxis.set_minor_locator((AutoMinorLocator(5)))

add_minor_ticks(scat_dist, 5)

scat_dist.plot(scat_bins[:-1], first_scat_dist, '-k', label=r'T$_1$')
scat_dist.plot(scat_bins[:-1], second_scat_dist, '--r', label=r'T$_2$')
scat_dist.plot(scat_bins[:-1], third_scat_dist, ':b', label=r'T$_3$')

# scat_dist.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
# scat_dist.tick_params(bottom=True, top=True, left=True, right=True)

# compton.plot(compton_angl_real, compton_eneg_real, 'ok', markerfacecolor=(.5,.5,.5,.5))

# compt_2d.imshow(np.log(np.rot90(h_compt)), cmap='hot')
X, Y = np.meshgrid(x_edges, y_edges)
im = compt_2d.pcolormesh(X, Y, np.log(h_compt), cmap='hot')
compt_2d.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True);
compt_2d.tick_params(bottom=True, top=True, left=True, right=True)
compt_2d.xaxis.set_minor_locator((AutoMinorLocator(5)))
compt_2d.yaxis.set_minor_locator((AutoMinorLocator(5)))

x_slice = 227
y_sl_c = 218
y_sl_w = 20
y_sl_l = y_sl_c - y_sl_w
y_sl_h = y_sl_c + y_sl_w

compt_2d.vlines(X[0,x_slice],Y[y_sl_l,0],Y[y_sl_h,0])

compt_slice.plot(Y[y_sl_l:y_sl_h,0], h_compt[x_slice,y_sl_l:y_sl_h], label=('slice at ' + str(X[0,x_slice]) + ' radians'))
print(X[0,x_slice])
compt_slice.set_xlabel('Energy (keV)')
compt_slice.set_ylabel('Counts')
compt_slice.set_ylim(bottom=0)

compt_slice.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True);
compt_slice.tick_params(bottom=True, top=True, left=True, right=True)
compt_slice.xaxis.set_minor_locator((AutoMinorLocator(5)))
compt_slice.yaxis.set_minor_locator((AutoMinorLocator(5)))

#ay.set_xscale('log')
# X[180:250],
font_scale = 16

ax.set_xlabel('Scatters before escaping')
ax.set_ylabel('% of all gammas')
# ax.set_title('Number of scatters before a gamma escapes,\narriving into 30 cm of water at 45 degrees')
ay.set_xlabel('Gamma Energy at Escape (keV)')
ay.set_ylabel('% of all gammas/keV')
# ay.set_title('Energy of escaping gammas,\narriving into 30 cm of water at 45 degrees')
aw.set_xlabel('highest energy scatter number')
aw.set_ylabel('% of all scatters')
# aw.set_title('Brightest scatters')
an.set_xlabel('energy deposited in highest energy scatter (keV)')
an.set_ylabel('% of all scatters with given energy/keV')
# an.set_title('Energy deposited by brightest scatter')

scat_dist.set_xlabel('distance into water (cm)')
scat_dist.set_ylabel('% of events per cm')

compt_2d.set_xlabel('scatter angle (radians)')
compt_2d.set_ylabel('deposited energy (keV)')

ax.set_ylim(bottom=0.)
ax.set_xlim(left=0.)
ay.set_ylim(bottom=0.)
# ay.set_ylim(top=6.8)
aw.set_ylim(bottom=0.)
an.set_ylim(bottom=0.)
ay.set_xlim(left=0.,right=511.)
aw.set_xlim(left=0.)
an.set_xlim(left=0.)
scat_dist.set_ylim(bottom=0., top=10.)
scat_dist.set_xlim(left=0., right=29.)

# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
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
ax.legend()
ay.legend()
# aw.legend()
# an.legend()
scat_dist.legend()
plt.tick_params(bottom=True, top=True, left=True, right=True)
# plt.tick_params(which = 'minor',bottom=False, top=False, left=False, right=False)
# plt.colorbar(im)
# plt.show()
version = "3"
subversion = "c"
fig_ax.savefig('containment_figures/scatters_before_escape_v' + version + subversion)
fig_ay.savefig('containment_figures/energy_at_escape_v' + version + subversion)
fig_aw.savefig('containment_figures/scatter_with_highest_energy_v' + version + subversion)
fig_an.savefig('containment_figures/energy_of_brightest_v' + version + subversion)
fig_scat_dist.savefig('containment_figures/scatter_dist_v' + version + subversion)
fig_compt_2d.savefig('containment_figures/compton_snake_v' + version + subversion)
fig_compt_slice.savefig('containment_figures/compton_slice_v' + version + subversion)


# x[x >= 0]
