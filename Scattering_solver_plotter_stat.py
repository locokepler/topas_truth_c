from operator import rshift
from readline import append_history_file
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MaxNLocator, MultipleLocator, AutoMinorLocator)

def plot_prettier(dpi=200, fontsize=15):
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

def add_minor_ticks(plot, ticks, bot=True, tp=True, lft=True, rght=True, xticks=None, yticks=None):
	plot.tick_params(which = 'minor', bottom=bot, top=tp, left=lft, right=rght)
	plot.tick_params(bottom=True, top=True, left=True, right=True)
	if (bot or tp):
		if (xticks != None):
			plot.xaxis.set_minor_locator(AutoMinorLocator(int(xticks)))
		else:
			plot.xaxis.set_minor_locator(AutoMinorLocator(int(ticks)))
	if (lft or rght):
		if (yticks != None):
			plot.yaxis.set_minor_locator(AutoMinorLocator(int(yticks)))
		else:
			plot.yaxis.set_minor_locator(AutoMinorLocator(int(ticks)))		


# sphere/sphere_1sigma.debug
# lor_data/dot_9_x_stat.debug
# full_data = np.loadtxt("sphere/sphere_2sigma.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test2.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test3.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test4.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test5.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test6.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test7.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test8.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test9.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test10.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test11.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test12.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test13.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test14.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test15.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test16.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test17.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/error_test18.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/fixed_n.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/sigma_test_1.0.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("/home/kepler/pet_simulations/sphere_dot/1mil_dot_scattered/resultant.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/scatter_eng_10keV_cutoff.debug", delimiter=',', skiprows=1)
full_data = np.loadtxt("lor_data/scatter_eng_20keV_trigger.debug", delimiter=',', skiprows=1)
# full_data = np.loadtxt("lor_data/scatter_eng_100keV_trigger.debug", delimiter=',', skiprows=1)
#no_in_pat = full_data[full_data[:,1] == 0]
transverse_miss = full_data[:,2] #no_in_pat[:,2]
long_miss = full_data[:,3]
# print(relevent[0:10])
alpha_scat = full_data[:,4]
alpha_scat = np.append(alpha_scat,full_data[:,9])
beta_scat = full_data[:,5]
beta_scat = np.append(beta_scat, full_data[:,10])
gamma_scat = full_data[:,6]
gamma_scat = np.append(gamma_scat, full_data[:,11])
delta_scat = full_data[:,7]
delta_scat = np.append(delta_scat, full_data[:,12])



# alpha_scat -= 1
# beta_scat  -= 1
# gamma_scat -= 1
# delta_scat -= 1

tran_miss_clean = transverse_miss[transverse_miss > 0]
long_miss_clean = long_miss[long_miss > 0]

given_bins = np.logspace(np.log10(0.0001), np.log10(10), 100)

# hist, bins = np.histogram(relevent, bins=given_bins, range=(0,100))

hist_trans, bins_trans = np.histogram(tran_miss_clean, bins=40, range=(0,20))
hist_long, bins_long = np.histogram(long_miss_clean, bins=40, range=(0, 20))
log_long, log_bins_long = np.histogram(long_miss_clean, bins=given_bins, range=(.0001,10))

alpha_hist, alpha_bins = np.histogram(alpha_scat, bins=9, range=(1,10))
beta_hist, beta_bins = np.histogram(beta_scat, bins=9, range=(1,10))
gamma_hist, gamma_bins = np.histogram(gamma_scat, bins=9, range=(1,10))
delta_hist, delta_bins = np.histogram(delta_scat, bins=9, range=(1,10))


# hist_original = 100. * hist_original
# hist_new = 100. * hist_new

fig_ax, ax = plt.subplots()
fig_ay, ay = plt.subplots()
fig_az, az = plt.subplots()
fig_az2, az2 = plt.subplots()
plt.subplots_adjust(bottom=.12)

ax.plot(bins_trans[:-1], hist_trans, 'Dk', label = 'LOR miss transverse')
ay.plot(bins_long[:-1], hist_long, 'Dk', label = 'LOR miss longitudinal')
az.plot(log_bins_long[:-1], log_long, 'Dk', label = 'LOR miss distance longitudinal')
az2.plot(log_bins_long[:-1], np.log(log_long), 'Dk', label = 'LOR miss distance for 2 gamma solution')


ax.set_ylim(bottom=0.)
ax.set_xlim(left=0.)
ay.set_ylim(bottom=0.)
ay.set_xlim(left=0.)
az.set_ylim(bottom=0.)

# plt.xlabel('Miss Distance (cm)')
# plt.ylabel('Number of Occurances')
# plt.title('Historgram of line of responce miss distance for scatters')
# plt.xlim(40, 160)
# plt.ylim(0, 0.03)
plt.grid(False)
from matplotlib.ticker import AutoMinorLocator
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('miss distance (cm)')
ay.set_xlabel('miss distance (cm)')
az.set_xscale('log')
az.set_xlabel('miss distance (cm)')
az2.set_xscale('log')
az2.set_yscale('log')
plt.tick_params(bottom=True, top=True, left=True, right=True)
plt.tick_params(which = 'minor',bottom=True, top=True, left=True, right=True)
ax.legend()
ay.legend()
az.legend()

fig_r_v_t, r_vs_t = plt.subplots()
plt.subplots_adjust(bottom=.12, left=0.2)

r_vs_t.set_ylim(bottom=0.)
r_vs_t.set_ylim(top=100.)
norm = np.shape(alpha_scat)[0] / 100.
r_vs_t.errorbar(alpha_bins[:-1], alpha_hist/norm, yerr=np.sqrt(alpha_hist)/norm, marker='D', linestyle='', label = r'T for R1')
r_vs_t.errorbar(beta_bins[:-1],  beta_hist/norm,  yerr=np.sqrt(beta_hist)/norm,  marker='s', linestyle='', label = r'T for R2')
r_vs_t.errorbar(gamma_bins[:-1], gamma_hist/norm, yerr=np.sqrt(gamma_hist)/norm, marker='o', linestyle='', label = r'T for R3')
r_vs_t.errorbar(delta_bins[:-1], delta_hist/norm, yerr=np.sqrt(delta_hist)/norm, marker='P', linestyle='', label = r'T for R4')
r_vs_t.legend()

r_vs_t.set_xlabel('True Scatter Number (T)')
r_vs_t.set_ylabel('% of Gammas')
r_vs_t.set_xlim(left=0.)

add_minor_ticks(r_vs_t, 5, xticks=2)

print("RvT values, rows are R, cols are T")
print(alpha_hist/norm)
print(beta_hist/norm)
print(gamma_hist/norm)
print(delta_hist/norm)

# information on second best guess deltas
# first_guess = full_data[:,4]
# second_guess = full_data[:,5]
# guesses = first_guess
# guesses = np.append(guesses,second_guess)
# guesses = guesses[guesses >= 0]

# log_guess_bins = np.logspace(np.log10(0.1), np.log10(100), 30)

# guess_hist, guess_bins = np.histogram(guesses, bins=log_guess_bins)#, range= (0,50))

# print(np.sum(guess_hist) / len(guesses))

# guess_hist = guess_hist / len(guesses)

# fig, guess_fig = plt.subplots()
# fig, scat_num = plt.subplots()

# guess_fig.plot(guess_bins[:-1], guess_hist, 'Dk', label = 'Difference between best and second best guess')
# guess_fig.set_ylim(bottom=0)

# scat_num.plot(alpha_bins[:-1], alpha_hist, 'sk', label='alpha')
# scat_num.plot(beta_bins[:-1], beta_hist, 'ok', label='beta')
# scat_num.set_xlim(left=0, right = 10)
# scat_num.set_ylim(bottom=0)
# scat_num.legend()

# # guess_fig.set_xlim(left=0)
# guess_fig.set_xscale('log')
# guess_fig.set_xlabel('energy delta (keV)')
# guess_fig.set_ylabel('fraction of solutions')
# guess_fig.legend()


fig_r_v_t.savefig('r_vs_t')
plt.show()
half_norm = 0.5 * norm
# print('Num of in patient scatters: ' + str(np.sum(full_data[:,1])/half_norm)) # number of in patient scatters

one_done = np.logical_or(full_data[:,4]!=-1, full_data[:,9]!=-1)
# print('Num of reconstructions of individual gammas: ' + str(np.sum(one_done)/half_norm))
# print('One good scatter in patient rate: ' + str(np.sum(one_done * full_data[:,1])/half_norm))

# percent of all events with both gammas sucessfully reconstructed
both_done = np.logical_and(full_data[:,4]!=-1,full_data[:,9]!=-1)
num_both_done = np.sum(both_done)
print('Num of reconstructions of both gammas: ' + str(num_both_done/half_norm))



# percent of all gammas that had a sucessful reconstruction
one_good = np.logical_or(full_data[:,4]==1, full_data[:,9]==1)
# print('Num of good reconstructions of individual gammas: ' + str(np.sum(one_good)/half_norm))
# print('One good scatter in patient rate: ' + str(np.sum(one_good * full_data[:,1])/half_norm))

# percent of all events with both gammas sucessfully reconstructed
good_a = np.logical_and(full_data[:,4]==1,full_data[:,9]==1)
print('Num of good reconstructions of both gammas: ' + str(np.sum(good_a)/half_norm))
# print('Two good scatter in patient rate: ' + str(np.sum(good_a * full_data[:,1])/half_norm))
print('Two scatter in patient rate: ' + str((100. * np.sum(both_done * full_data[:,1]))/num_both_done))

#snr: both good / total reconstructions
print('SNR (good/(done - good)): ' + str(np.sum(np.logical_and(good_a, full_data[:,1]==0)) / (num_both_done - np.sum(good_a))))

# x[x >= 0]
