# A script to plot the tracer density from the 
# Bryan Fritsch bubble test

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

set_tomplot_style()

order = 1
dxz = 50

# For results on the previous branch, which used Td with interpolate
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/bryan_fritsch'
#results_dir = '/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/bryan_fritsch/'

# Results in another previous branch
#results_dir = '/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/bryan_fritsch/'
#base_file_name = 'bryan_fritsch_linear_mX0_16x1'

# Results in the current branch
results_dir = ''
base_file_name = 'bryan_fritsch_linear_mX0_16x1'
#base_file_name = 'bryan_fritsch_feb18_2026'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/bryan_fritsch_order_{order}_tracer_density.png'

# Advective result:
adv_dirname = f'{base_file_name}_advective_order_{order}dxz{dxz}'
nc = Dataset(f'{results_dir}results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])

T_wv = nc.groups['TracerDensity_water_vapour_rho']
Twv = T_wv['total'][:]
T_cw = nc.groups['TracerDensity_cloud_water_rho']
Tcw = T_cw['total'][:]

Td_adv = Twv + Tcw
Td_adv_diff = np.abs(Td_adv - Td_adv[0])/Td_adv[0]

# Conservative tracer result:
con_dirname = f'{base_file_name}_conservative_order_{order}dxz{dxz}'
conservative_nc = Dataset(f'{results_dir}results/{con_dirname}/diagnostics.nc')

T_wv = conservative_nc.groups['TracerDensity_water_vapour_rho']
Twv = T_wv['total'][:]
T_cw = conservative_nc.groups['TracerDensity_cloud_water_rho']
Tcw = T_cw['total'][:]

Td_con = Twv + Tcw
Td_con_diff = np.abs(Td_con - Td_con[0])/Td_con[0]

print('End Td diff for advective', Td_adv_diff[-1])
print('End Td diff for conservative', Td_con_diff[-1])

# Plot these tracer densities
fig, ax = plt.subplots(1, 1, figsize=(7, 5))
plt.semilogy(times, Td_adv_diff, label='advective', c='r')
plt.semilogy(times, Td_con_diff, label='conservative', c='b')
plt.xlabel('Time (s)', size = 16)
plt.ylabel('Relative change in tracer density', size=16)
plt.gca().ticklabel_format(axis='x',useMathText=True)
plt.legend(loc='lower center', prop={'size': 16}, bbox_to_anchor=(0.5, -0.4))
ax.set_xlim([0,1000])
ax.set_ylim([1e-16, 1e-5])

#tomplot_legend_ax(ax, location='bottom')

print(f'Saving figure to {plot_name}')
plt.savefig(plot_name, bbox_inches='tight')
plt.close()