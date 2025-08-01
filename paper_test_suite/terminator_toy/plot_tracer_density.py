# A script to plot the tracer density from the 
# Terminator Toy chemsitry test

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

order = 1
dxz = 50

ncells_per_edge = 24

quantity = 'tracer_density'

# Specific names for the advective and conservative results:
adv_extra_name = 'analyt_ref_Td_solve_'
con_extra_name = 'aug1_zerocrops_nochemlim2'
save_name = 'aug1_zerocrops_nochemlim2'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/terminator_toy_{save_name}_{quantity}_over_time_{dxz}.png'

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

# Conservative result:
# NEED TO CHANGE BACK TO CONERSVATIVE!!
con_dirname = f'terminator_toy_{con_extra_name}conservative_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{con_dirname}/diagnostics.nc')

# Conservative
Td_con_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_con_X = Td_con_X['total'][:]
Td_con_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_con_X2 = Td_con_X2['total'][:]
Td_con = Td_con_X + 2*Td_con_X2

Td_con_diff = np.abs(Td_con - Td_con[0])/Td_con[0]

# Plot these tracer densities
plt.figure()
plt.plot(times, Td_adv_diff, label='advective')
plt.plot(times, Td_con_diff, label='conservative')
plt.xlabel('Time (s)')
plt.ylabel('$|T_{d}(t) - T_{d}(0)| / T_{d}(0)$')
plt.title(f'Species conservation in the Terminator Toy')
plt.legend()

#print(f'Saving figure to {plot_name}')
#plt.savefig(plot_name)
#plt.close()

#Check both have the same start densities:
print('Diff in initial Td between results is ', np.abs(Td_con[0]-Td_adv[0]))

print('End Td diff for advective', Td_adv_diff[-1])
print('End Td diff for conservative', Td_con_diff[-1])


#Log graph
plt.figure()
#plt.semilogy(times, Td_adv_diff1, label='advective')
plt.plot(times, Td_adv_diff, label='advective')
plt.semilogy(times, Td_con_diff, label='conservative')
plt.xlabel('Time (s)')
plt.ylabel('$|T_{d}(t) - T_{d}(0)| / T_{d}(0)$')
plt.title(f'Species conservation in the Terminator Toy')
plt.legend()

plot_name_semilogy = f'{plot_dir}/terminator_toy_{save_name}_Td_over_time_dxz{dxz}_semilogy.png'

print(f'Saving figure to {plot_name_semilogy}')
plt.savefig(plot_name_semilogy)
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

print(f'Saving figure to {plot_name_min}')
plt.savefig(plot_name_min)
plt.close()