import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

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
			

# source_file = "lor_data/LAB_sensitivity_1mm_1keV_500ps"
# source_file = "lor_data/LAB_sensitivity_0.1mm_1keV_500ps_split"
# source_file = "lor_data/LAB_sensitivity_0.32mm_1keV_500ps_split"
# source_file = "lor_data/LAB_sensitivity_1mm_100keV_500ps"
# source_file = "lor_data/LAB_sensitivity_1mm_10keV_500ps"
# source_file = "lor_data/LAB_sensitivity_3.2mm_1keV_500ps_split"
# source_file = "lor_data/LAB_sensitivity_1mm_5keV_500ps"
# source_file = "lor_data/LAB_sensitivity_1mm_10keV_500ps"
# source_file = "lor_data/LAB_sensitivity_1mm_50keV_500ps"
# source_file = "lor_data/LAB_sensitivity_1mm_100keV_500ps"
# source_file = "lor_data/LAB_sensitivity_1mm_1keV_50ps"
# source_file = "lor_data/LAB_sensitivity_1mm_1keV_100ps"
# source_file = "lor_data/LAB_sensitivity_1mm_1keV_1000ps"
# source_file = "lor_data/LAB_sensitivity_1mm_1keV_5000ps"

# source_file = "lor_data/70cm_sensitivity_test"
source_file = "lor_data/LYSO_sensitivity_energy_10"
# source_lors = np.loadtxt(source_file + '.misID', delimiter=',')
source_lors = np.loadtxt(source_file + '.lor', delimiter=',')

# source_lors[:,0] hist_nums 
# source_lors[:,1:3] centers
# source_lors[:,4:6] vectors

cutoff_cm = 0.2
print(np.shape(source_lors))

# find how close the LOR goes to 0,0,0
# this can be done by center position cross unit vector of direction

unit_velocity = np.array([source_lors[:,4], source_lors[:,5], source_lors[:,6]])/np.hypot(np.hypot(source_lors[:,4],source_lors[:,5]),source_lors[:,6])
print(np.shape(unit_velocity))
center = np.array([source_lors[:,1], source_lors[:,2], source_lors[:,3]])
print(np.shape(center))
crosses = np.cross(unit_velocity,center, axis=0)
distance = np.hypot(np.hypot(crosses[0,:], crosses[1,:]), crosses[2,:])
print(np.shape(distance))

log_bins = np.logspace(-3,2, 5 * 10)
dist_hist, dist_bins = np.histogram(distance, bins=log_bins)
# print(dist_bins)

fig_ax, ax = plt.subplots()
plt.subplots_adjust(bottom=0.13,top=0.88,left=0.21,right=0.96)

ax.bar(dist_bins[:-1],dist_hist/np.sum(dist_hist), width=np.diff(dist_bins))
add_minor_ticks(ax, 5)
ax.set_xlabel("miss distance (cm)")
ax.set_ylabel("fraction")
ax.set_xscale('log')
ax.vlines(cutoff_cm, 0, np.amax(dist_hist)/np.sum(dist_hist), color='green', linestyle='dashed')

fig_ax.savefig('lor_miss_split.png')

good_lors = source_lors[distance < cutoff_cm]
bad_lors = source_lors[distance >= cutoff_cm]
print(np.shape(good_lors))
print(np.shape(bad_lors))

print(np.shape(good_lors)[0] / np.shape(source_lors)[0])

np.savetxt(source_file + ".good", good_lors, fmt='%d,%e,%e,%e,%e,%e,%e,%e,%e')
np.savetxt(source_file + ".bad", bad_lors, fmt='%d,%e,%e,%e,%e,%e,%e,%e,%e')
