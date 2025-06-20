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
adv_ref_name = 'analyt_ref_Td_solve'
adv_extra_name = 'mmr_twolambda'
con_ref_name = 'ref_analyt_forced'
con_extra_name = 'mmr_twolambda'
save_name = 'two_lambda'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/terminator_toy_{save_name}_with_ref_Td_over_time_{dxz}.png'

# Advective ref result:
adv_ref_dirname = f'terminator_toy_{adv_ref_name}_advective_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{adv_ref_dirname}/diagnostics.nc')

# Advective ref
Td_adv_ref_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_adv_ref_X = Td_adv_ref_X['total'][:]
Td_adv_ref_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_adv_ref_X2 = Td_adv_ref_X2['total'][:]
Td_adv_ref = Td_adv_ref_X + 2*Td_adv_ref_X2
Td_adv_ref_diff = np.abs(Td_adv_ref - Td_adv_ref[0])/Td_adv_ref[0]

# Advective MMR result:
adv_dirname = f'terminator_toy_{adv_extra_name}_advective_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{adv_dirname}/diagnostics.nc')

# Advective MMR
Td_adv_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_adv_X = Td_adv_X['total'][:]
Td_adv_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_adv_X2 = Td_adv_X2['total'][:]
Td_adv = Td_adv_X + 2*Td_adv_X2
Td_adv_diff = np.abs(Td_adv - Td_adv[0])/Td_adv[0]

################################
# Conservative ref result:
con_ref_dirname = f'terminator_toy_{con_ref_name}_conservative_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{con_ref_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])
print(nc.groups)

# Conservative ref
Td_con_ref_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_con_ref_X = Td_con_ref_X['total'][:]
Td_con_ref_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_con_ref_X2 = Td_con_ref_X2['total'][:]
Td_con_ref = Td_con_ref_X + 2*Td_con_ref_X2
Td_con_ref_diff = np.abs(Td_con_ref - Td_con_ref[0])/Td_con_ref[0]

# Conservative MMR result:
con_dirname = f'terminator_toy_{con_extra_name}_conservative_ncells_{ncells_per_edge}'
nc = Dataset(f'results/{con_dirname}/diagnostics.nc')

# Conservative MMR
Td_con_X = nc.groups['TracerDensity_X_tracer_rho_d']
Td_con_X = Td_con_X['total'][:]
Td_con_X2 = nc.groups['TracerDensity_X2_tracer_rho_d']
Td_con_X2 = Td_con_X2['total'][:]
Td_con = Td_con_X + 2*Td_con_X2

Td_con_diff = np.abs(Td_con - Td_con[0])/Td_con[0]

# Plot these tracer densities
plt.figure()
plt.plot(times, Td_adv_ref_diff, label='advective reference')
plt.plot(times, Td_adv_diff, label='advective MMR ')
plt.plot(times, Td_con_ref_diff, label='conservative reference')
plt.plot(times, Td_con_diff, label='conservative MMR')
plt.xlabel('Time (s)')
plt.ylabel('$|T_{d}(t) - T_{d}(0)| / T_{d}(0)$')
plt.title(f'Species conservation in the Terminator Toy')
plt.legend()

print(f'Saving figure to {plot_name}')
plt.savefig(plot_name)
plt.close()

#Check both have the same start densities:
print('Diff in initial Td between results is ', np.abs(Td_con[0]-Td_adv[0]))

print('End Td diff for advective reference', Td_adv_ref_diff[-1])
print('End Td diff for advective MMR', Td_adv_diff[-1])
print('End Td diff for conservative reference', Td_con_ref_diff[-1])
print('End Td diff for conservative MMR', Td_con_diff[-1])


#Log graph
plt.figure()
plt.semilogy(times, Td_adv_ref_diff, label='advective reference')
plt.semilogy(times, Td_adv_diff, label='advective MMR')
plt.semilogy(times, Td_con_ref_diff, label='conservative reference')
plt.semilogy(times, Td_con_diff, label='conservative MMR')
plt.xlabel('Time (s)')
plt.ylabel('$|T_{d}(t) - T_{d}(0)| / T_{d}(0)$')
plt.title(f'Species conservation in the Terminator Toy')
plt.legend()

plot_name_semilogy = f'{plot_dir}/terminator_toy_{save_name}_with_ref_Td_over_time_dxz{dxz}_semilogy.png'

print(f'Saving figure to {plot_name_semilogy}')
plt.savefig(plot_name_semilogy)
plt.close()