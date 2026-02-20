'''
Plot the tracer density over time from the simple
for the MMR limiter
'''

import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import matplotlib
import matplotlib.colors as colors
import cartopy.crs as ccrs
import os
from os.path import abspath, dirname
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_ax, tomplot_field_title, tomplot_contours,
    extract_gusto_coords, extract_gusto_field, reshape_gusto_data
)

########################################

set_tomplot_style()

no_lim_nc = Dataset(f'{abspath(dirname(__file__))}/results/mmr_test_no_limiter_sphere_cylinder/diagnostics.nc')
mmr_lim_nc = Dataset(f'{abspath(dirname(__file__))}/results/mmr_test_mmr_limiter_sphere_cylinder/diagnostics.nc')

extra_name = '_cylinder'

figure_stem = f'{abspath(dirname(__file__))}/figures/'

no_lim_td = no_lim_nc.groups['TracerDensity_m_X_rho_d']
mmr_lim_td = mmr_lim_nc.groups['TracerDensity_m_X_rho_d']

no_lim_mx_error = no_lim_nc.groups['m_X_error']['l2'][-1]
mmr_lim_mx_error = mmr_lim_nc.groups['m_X_error']['l2'][-1]

print('No limiter m_X error is', no_lim_mx_error)
print('MMR limiter m_X error is', mmr_lim_mx_error)

no_lim_td = no_lim_td['total'][:]
mmr_lim_td = mmr_lim_td['total'][:]

print('Initial tracer density is', no_lim_td[0])

no_lim_td_diff = np.abs(no_lim_td - no_lim_td[0])/no_lim_td[0]
mmr_lim_td_diff = np.abs(mmr_lim_td - mmr_lim_td[0])/mmr_lim_td[0]

times = np.asarray(no_lim_nc['time'])
time_days = times/60/60/24

fig, ax = plt.subplots(1,1,figsize=(7,5))
ax.semilogy(time_days, no_lim_td_diff, label='No limiter')
ax.semilogy(time_days, mmr_lim_td_diff, label='MMR limiter')
ax.set_ylabel('Relative change in tracer density',size=16)
ax.set_xlabel('Time (days)', size=16)
plt.legend(loc='lower center', prop={'size': 16}, bbox_to_anchor=(0.5, -0.4))
ax.set_xlim([0,12])
ax.set_ylim([1e-16,1e-11])
plt.grid()

savename = f'{figure_stem}mmr_test_sphere_tracer_density{extra_name}.jpg'
plt.savefig(savename, bbox_inches='tight')
print(f'saved figure to {savename}')

plt.close()