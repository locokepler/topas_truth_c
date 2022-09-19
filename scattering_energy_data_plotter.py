import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MaxNLocator, MultipleLocator, AutoMinorLocator)


# full_data = np.loadtxt("sphere/sphere_2sigma_skp2.eng", delimiter=',', skiprows=0)
full_data = np.loadtxt("lor_data/fixed_n.eng", delimiter=',', skiprows=0)

plot_version = 3

def add_minor_ticks(plot, ticks):
	plot.tick_params(which = 'minor', bottom=True, top=True, left=True, right=True)
	plot.tick_params(bottom=True, top=True, left=True, right=True)
	plot.xaxis.set_minor_locator((AutoMinorLocator(int(ticks))))
	plot.yaxis.set_minor_locator((AutoMinorLocator(int(ticks))))


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

denomenator = np.shape(full_data)[0] * 2 # number of gammas as the data source

alpha_eng = full_data[:,0]
alpha_eng = np.append(alpha_eng, full_data[:,5])
beta_eng  = full_data[:,1]
beta_eng  = np.append(beta_eng, full_data[:,6])
gamma_eng = full_data[:,2]
gamma_eng = np.append(gamma_eng, full_data[:,7])
delta_eng = full_data[:,3]
delta_eng = np.append(delta_eng, full_data[:,8])

clean_alpha = alpha_eng[alpha_eng >= 0]
clean_beta  = beta_eng[beta_eng >= 0]
clean_gamma = gamma_eng[gamma_eng >= 0]
clean_delta = delta_eng[delta_eng >= 0]

energy_bins = 20
energy_range = (0, 340)

norm = 100 * (energy_bins / (denomenator * (energy_range[1] - energy_range[0])))
# normalization of histogram per unit (keV)

alpha_hist, alpha_bins = np.histogram(clean_alpha, bins=energy_bins, range=energy_range)
beta_hist,  beta_bins  = np.histogram(clean_beta,  bins=energy_bins, range=energy_range)
gamma_hist, gamma_bins = np.histogram(clean_gamma, bins=energy_bins, range=energy_range)
delta_hist, delta_bins = np.histogram(clean_delta, bins=energy_bins, range=energy_range)

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=.15)

# ax.plot(alpha_bins[:-1], alpha_hist * norm, '-', label = r'E$_{deposit}$ for R1')
ax.errorbar(alpha_bins[:-1], alpha_hist * norm, yerr=np.sqrt(alpha_hist) * norm, marker='D', label = r'E$_{deposit}$ for R$_1$')
# ax.plot(beta_bins[:-1],  beta_hist * norm,  '--', label = r'E$_{deposit}$ for R2')
ax.errorbar(beta_bins[:-1], beta_hist * norm, yerr=np.sqrt(beta_hist) * norm, marker='s', label = r'E$_{deposit}$ for R$_2$')
# ax.plot(gamma_bins[:-1], gamma_hist * norm, '-.', label = r'E$_{deposit}$ for R3')
ax.errorbar(gamma_bins[:-1], gamma_hist * norm, yerr=np.sqrt(gamma_hist) * norm, marker='o', label = r'E$_{deposit}$ for R$_3$')
# ax.plot(delta_bins[:-1], delta_hist * norm, ':', label = r'E$_{deposit}$ for R4')
ax.errorbar(delta_bins[:-1], delta_hist * norm, yerr=np.sqrt(delta_hist) * norm, marker='p', label = r'E$_{deposit}$ for R$_4$')

# ax.errorbar(alpha_bins[:-1], np.log(alpha_hist * norm), yerr=np.log(np.sqrt(alpha_hist) * norm), marker='D', label = r'E$_{deposit}$ for R1')
# ax.errorbar(beta_bins[:-1], np.log(beta_hist * norm), yerr=np.log(np.sqrt(beta_hist) * norm), marker='s', label = r'E$_{deposit}$ for R2')
# ax.errorbar(gamma_bins[:-1], np.log(gamma_hist * norm), yerr=np.log(np.sqrt(gamma_hist) * norm), marker='o', label = r'E$_{deposit}$ for R3')
# ax.errorbar(delta_bins[:-1], np.log(delta_hist * norm), yerr=np.log(np.sqrt(delta_hist) * norm), marker='p', label = r'E$_{deposit}$ for R4')


ax.set_ylim(bottom=0.)
ax.set_xlim(left=0.)

add_minor_ticks(ax, 5)


ax.set_xlabel(r'E$_{deposit}$ (keV)')
ax.set_ylabel('% of gammas/keV')

ax.legend()

plt.savefig('R_eng_errorbar_v'+str(plot_version))

n1_eng = full_data[:,10]
n1_eng = np.append(n1_eng, full_data[:,15])
n2_eng  = full_data[:,11]
n2_eng  = np.append(n2_eng, full_data[:,16])
n3_eng = full_data[:,12]
n3_eng = np.append(n3_eng, full_data[:,17])
n4_eng = full_data[:,13]
n4_eng = np.append(n4_eng, full_data[:,18])

clean_1n = n1_eng[n1_eng >= 0]
clean_2n = n2_eng[n2_eng >= 0]
clean_3n = n3_eng[n3_eng >= 0]
clean_4n = n4_eng[n4_eng >= 0]

n1_hist, n1_bins = np.histogram(clean_1n, bins=energy_bins, range=energy_range)
n2_hist, n2_bins = np.histogram(clean_2n, bins=energy_bins, range=energy_range)
n3_hist, n3_bins = np.histogram(clean_3n, bins=energy_bins, range=energy_range)
n4_hist, n4_bins = np.histogram(clean_4n, bins=energy_bins, range=energy_range)

fig, ay = plt.subplots()
plt.subplots_adjust(bottom=.15)


# ay.plot(n1_bins[:-1], n1_hist * norm, '-', label = r'E$_{deposit}$ for T1')
# ay.plot(n2_bins[:-1], n2_hist * norm,  '--', label = r'E$_{deposit}$ for T2')
# ay.plot(n3_bins[:-1], n3_hist * norm, '-.', label = r'E$_{deposit}$ for T3')
# ay.plot(n4_bins[:-1], n4_hist * norm, ':', label = r'E$_{deposit}$ for T4')

ay.errorbar(n1_bins[:-1], n1_hist * norm, yerr=np.sqrt(n1_hist) * norm, marker='D', label = r'E$_{deposit}$ for T1')
ay.errorbar(n2_bins[:-1], n2_hist * norm, yerr=np.sqrt(n2_hist) * norm, marker='s', label = r'E$_{deposit}$ for T2')
ay.errorbar(n3_bins[:-1], n3_hist * norm, yerr=np.sqrt(n3_hist) * norm, marker='o', label = r'E$_{deposit}$ for T3')
ay.errorbar(n4_bins[:-1], n4_hist * norm, yerr=np.sqrt(n4_hist) * norm, marker='p', label = r'E$_{deposit}$ for T4')


add_minor_ticks(ay, 5)

ay.set_ylim(bottom=0.)
ay.set_xlim(left=0.)


ay.set_xlabel(r'E$_{deposit}$ (keV)')
ay.set_ylabel('% of gammas/keV')

ay.legend()


plt.savefig('T_eng_errorbar_v'+str(plot_version))
#plt.show()
