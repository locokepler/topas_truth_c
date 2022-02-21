import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numeric import full

full_data = np.loadtxt("test", delimiter=',', skiprows=1)
#no_in_pat = full_data[full_data[:,1] == 0]
original_algo = full_data[:,2] #no_in_pat[:,2]
new_algo = full_data[:,3]
# print(relevent[0:10])
alpha_scat = full_data[:,6]
alpha_scat = np.append(alpha_scat,full_data[:,8])
beta_scat = full_data[:,7]
beta_scat = np.append(beta_scat, full_data[:,9])

a_orginal = original_algo[original_algo > 0]
a_new = new_algo[new_algo > 0]

given_bins = np.logspace(np.log10(0.0001), np.log10(10), 100)

# hist, bins = np.histogram(relevent, bins=given_bins, range=(0,100))

hist_original, bins_original = np.histogram(a_orginal, bins=40, range=(0,80))
hist_new, bins_new = np.histogram(a_new, bins=40, range=(0, 80))
log_new, log_bins_new = np.histogram(a_new, bins=given_bins, range=(.0001,10))

alpha_hist, alpha_bins = np.histogram(alpha_scat, bins=10, range=(1,10))
beta_hist, beta_bins = np.histogram(beta_scat, bins=10, range=(1,10))

# hist_original = 100. * hist_original
# hist_new = 100. * hist_new

fig, ax = plt.subplots()
fig, ay = plt.subplots()
fig, az = plt.subplots()
fig, az2 = plt.subplots()

ax.plot(bins_original[:-1], hist_original, 'Dk', label = 'LOR miss for 1 gamma solution')
ay.plot(bins_new[:-1], hist_new, 'Dk', label = 'LOR miss for 2 gamma solution')
az.plot(log_bins_new[:-1], log_new, 'Dk', label = 'LOR miss distance for 2 gamma solution')
az2.plot(log_bins_new[:-1], np.log(log_new), 'Dk', label = 'LOR miss distance for 2 gamma solution')


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

# information on second best guess deltas
first_guess = full_data[:,4]
second_guess = full_data[:,5]
guesses = first_guess
guesses = np.append(guesses,second_guess)
guesses = guesses[guesses >= 0]

log_guess_bins = np.logspace(np.log10(0.1), np.log10(100), 30)

guess_hist, guess_bins = np.histogram(guesses, bins=log_guess_bins)#, range= (0,50))

print(np.sum(guess_hist) / len(guesses))

guess_hist = guess_hist / len(guesses)

fig, guess_fig = plt.subplots()
fig, scat_num = plt.subplots()

guess_fig.plot(guess_bins[:-1], guess_hist, 'Dk', label = 'Difference between best and second best guess')
guess_fig.set_ylim(bottom=0)

scat_num.plot(alpha_bins[:-1], alpha_hist, 'sk', label='alpha')
scat_num.plot(beta_bins[:-1], beta_hist, 'ok', label='beta')
scat_num.set_xlim(left=0, right = 10)
scat_num.set_ylim(bottom=0)
scat_num.legend()

# guess_fig.set_xlim(left=0)
guess_fig.set_xscale('log')
guess_fig.set_xlabel('energy delta (keV)')
guess_fig.set_ylabel('fraction of solutions')
guess_fig.legend()


plt.show()



# x[x >= 0]
