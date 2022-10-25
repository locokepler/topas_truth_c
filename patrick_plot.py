from fileinput import filename
from sre_parse import fix_flags
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

# file_name = "lor_data/in_pat_1mil_0_8_no_time"
# file_name = "/home/kepler/pet_simulations/sphere_dot/1mil_dot_scattered/resultant"
# file_name = "lor_data/types_no_IPS1"
# file_name = "lor_data/types_no_IPS1_3"
# file_name = "lor_data/types_no_IPS_OOPW_0"
# file_name = "lor_data/types_IPS_OOPW_0"
file_name = "lor_data/types_no_IPS_keep_all"
# file_name = "lor_data/types_IPS"
print(file_name)

gold_data = np.loadtxt(file_name + ".gold", delimiter=',')
silver_data = np.loadtxt(file_name + ".silver", delimiter=',')
lead_data = np.loadtxt(file_name + ".lead", delimiter=',')

# reduce down to just the x axis center value
gold_x = gold_data[:,1]
silver_x = silver_data[:,1]
lead_x = lead_data[:,1]

num_of_bins = 81

gold_hist, gold_bins = np.histogram(gold_x, bins=num_of_bins, range=(-10,10))
silver_hist, silver_bins = np.histogram(silver_x, bins=num_of_bins, range=(-10,10))
lead_hist, lead_bins = np.histogram(lead_x, bins=num_of_bins, range=(-10,10))

fig_x, x_hist = plt.subplots()
plt.subplots_adjust(bottom=0.12, left=0.18, right=0.92)

bin_offset = (gold_bins[1] - gold_bins[0]) / 2

x_hist.set_yscale('log')
x_hist.plot(gold_bins[:-1] + bin_offset,gold_hist, label='gold events')
x_hist.plot(gold_bins[:-1] + bin_offset,silver_hist, label='silver events')
x_hist.plot(gold_bins[:-1] + bin_offset,lead_hist, label='lead events')

x_hist.set_xlabel("x distance from center")
x_hist.set_ylabel("number of LORs")

x_hist.tick_params(which = 'minor', bottom=True, top=True)
x_hist.tick_params(which = 'major', bottom=True, top=True)
x_hist.xaxis.set_minor_locator(AutoMinorLocator(5))

x_hist.legend()

fig_x_lin, x_hist_lin = plt.subplots()
plt.subplots_adjust(bottom = 0.12, left=0.18, right=0.92)
x_hist_lin.plot(gold_bins[:-1] + bin_offset,gold_hist, label='gold events')
x_hist_lin.plot(gold_bins[:-1] + bin_offset,silver_hist, label='silver events')
x_hist_lin.plot(gold_bins[:-1] + bin_offset,lead_hist, label='lead events')
x_hist_lin.set_ylim(bottom=0.0)

add_minor_ticks(x_hist_lin,5)
x_hist_lin.legend()

print("Gold lors: " + str(len(gold_x)))
print("Silver lors: " + str(len(silver_x)))
print("Lead lors: " + str(len(lead_x)))

fig_x.savefig('patrick_plot_v2')
fig_x_lin.savefig('patrick_plot_lin_v2')