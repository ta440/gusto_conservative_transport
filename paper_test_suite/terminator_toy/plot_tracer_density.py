# A script to plot the tracer density from the 
# Terminator Toy chemsitry test

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

set_tomplot_style()

ncells_per_edge = 24

quantity = 'tracer_density'

# Specific names for the advective and conservative results:
#adv_extra_name = 'analyt_ref_Td_solve_'
#con_extra_name = 'aug1_zerocrops_nochemlim2'
#save_name = 'aug1_zerocrops_nochemlim2'

#adv_extra_name = 'dg1_lim_'
#con_extra_name = ''
#save_name = ''

#adv_extra_name = 'no_phys_no_lim_'
#con_extra_name = 'no_phys_'
#save_name = 'no_phys_'

adv_extra_name = 'dg1_lim_'
con_extra_name = 'feb3_'
save_name = 'feb3_'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/terminator_toy_{save_name}{quantity}_over_time_ncells_{ncells_per_edge}.png'

# Advective result:
adv_dirname = f'terminator_toy_{adv_extra_name}advective_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])
print(nc.groups)

# Advective
Td_adv_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_adv_X = Td_adv_X['total'][:]
Td_adv_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_adv_X2 = Td_adv_X2['total'][:]
Td_adv = Td_adv_X + 2*Td_adv_X2

Td_adv_diff = np.abs(Td_adv - Td_adv[0])/Td_adv[0]
Td_adv_X_diff = np.abs(Td_adv_X - Td_adv_X[0])/Td_adv_X[0]
Td_adv_X2_diff = np.abs(Td_adv_X2 - Td_adv_X2[0])/Td_adv_X2[0]

# Mixing ratios:
X_adv = nc.groups['X_tracer']['total'][:]
X2_adv = nc.groups['X2_tracer']['total'][:]

X_adv_diff = X_adv - X_adv[0]
X2_adv_diff = X2_adv - X2_adv[0]

print('End X for advective scheme', X_adv_diff[-1]/X_adv[0])
print('End X2 for advective scheme', X2_adv_diff[-1]/X2_adv[0])

# Conservative result:
# NEED TO CHANGE BACK TO CONERSVATIVE!!
con_dirname = f'terminator_toy_{con_extra_name}conservative_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{con_dirname}/diagnostics.nc')

# Conservative
# Tracer densities
Td_con_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_con_X = Td_con_X['total'][:]
Td_con_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_con_X2 = Td_con_X2['total'][:]
Td_con = Td_con_X + 2*Td_con_X2

# Mixing ratios:
X_con = nc.groups['X_tracer']['total'][:]
X2_con = nc.groups['X2_tracer']['total'][:]

X_con_diff = X_con - X_con[0]
X2_con_diff = X2_con - X2_con[0]

print('End X for conservative scheme', X_con_diff[-1]/X_con[0])
print('End X2 for conservative scheme', X2_con_diff[-1]/X2_con[0])

Td_con_diff = np.abs(Td_con - Td_con[0])/Td_con[0]
Td_con_X_diff = np.abs(Td_con_X - Td_con_X[0])/Td_con_X[0]
Td_con_X2_diff = np.abs(Td_con_X2 - Td_con_X2[0])/Td_con_X2[0]

#Check both have the same start densities:
print('Diff in initial Td between results is ', np.abs(Td_con[0]-Td_adv[0]))

print('End Td diff for advective', Td_adv_diff[-1])
print('End Td diff for conservative', Td_con_diff[-1])

print('End Td X diff for advective', Td_adv_X_diff[-1])
print('End Td X diff for conservative', Td_con_X_diff[-1])

print('End Td X2 diff for advective', Td_adv_X2_diff[-1])
print('End Td X2 diff for conservative', Td_con_X2_diff[-1])

print(f'Advective dir: {adv_dirname}, {adv_dirname}')
print(f'Conservative dir: {con_dirname}, {con_dirname}')

times = np.asarray(nc['time'])
time_days = times/60/60/24


#Log graph
fig, ax = plt.subplots(1,1,figsize=(7,5))
ax.semilogy(time_days, Td_adv_diff, label='advective', c='r')
ax.semilogy(time_days, Td_con_diff, label='conservative', c='b')
#ax.semilogy(time_days, Td_adv_X_diff, label='advective X', linestyle='dashed', c='k')
#ax.semilogy(time_days, Td_con_X_diff, label='tracer conservative X', linestyle='dashed', c='k')
#ax.semilogy(time_days, Td_adv_X2_diff, label='advective X2',linestyle='dotted')
#ax.semilogy(time_days, Td_con_X2_diff, label='tracer conservative X2',linestyle='dotted')
ax.set_xlabel('Time (days)', size=16)
ax.set_ylabel('Relative change in tracer density', size=16)
plt.gca().ticklabel_format(axis='x',useMathText=True)
plt.legend(loc='lower center', prop={'size': 16}, bbox_to_anchor=(0.5, -0.4))
ax.set_xlim([0,12])
ax.set_ylim([1e-16,1e-6])
plt.grid()


plot_name = f'{plot_dir}/terminator_toy_{save_name}Td_over_time_ncells_{ncells_per_edge}.png'

print(f'Saving figure to {plot_name}')
plt.savefig(plot_name, bbox_inches='tight')
plt.close()

# Bonus, plot the minimum value of the mixing ratios:
X_min = nc.groups['X_tracer']
X_min = X_min['min'][:]

X2_min = nc.groups['X2_tracer']
X2_min = X2_min['min'][:]

plt.figure()
plt.plot(times, X_min, label='X')
plt.plot(times, X2_min, label='X2')
plt.legend()

print('Min min value for X:', min(X_min))
print('Min min value for X2:', min(X2_min))

plot_name_min = f'{plot_dir}/terminator_toy_{save_name}_minvals.png'

#print(f'Saving figure to {plot_name_min}')
#plt.savefig(plot_name_min)
#plt.close()