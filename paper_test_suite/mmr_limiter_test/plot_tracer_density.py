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

no_lim_nc = Dataset(f'{abspath(dirname(__file__))}/results/mmr_test_no_limiter/diagnostics.nc')
dg1_lim_nc = Dataset(f'{abspath(dirname(__file__))}/results/mmr_test_standard_limiter/diagnostics.nc')
mmr_lim_nc = Dataset(f'{abspath(dirname(__file__))}/results/mmr_test_mmr_limiter5/diagnostics.nc')

figure_stem = f'{abspath(dirname(__file__))}/figures/'

no_lim_td = no_lim_nc.groups['TracerDensity_m_X_rho_d']
dg1_lim_td = dg1_lim_nc.groups['TracerDensity_m_X_rho_d']
mmr_lim_td = mmr_lim_nc.groups['TracerDensity_m_X_rho_d']

no_lim_mx_error = no_lim_nc.groups['m_X_error']['l2'][-1]
dg1_lim_mx_error = dg1_lim_nc.groups['m_X_error']['l2'][-1]
mmr_lim_mx_error = mmr_lim_nc.groups['m_X_error']['l2'][-1]

print('No limiter m_X error is', no_lim_mx_error)
print('Standard DG1 limiter m_X error is', dg1_lim_mx_error)
print('MMR limiter m_X error is', mmr_lim_mx_error)

no_lim_td = no_lim_td['total'][:]
dg1_lim_td = dg1_lim_td['total'][:]
mmr_lim_td = mmr_lim_td['total'][:]

no_lim_td_diff = np.abs(no_lim_td - no_lim_td[0])/no_lim_td[0]
dg1_lim_td_diff = np.abs(dg1_lim_td - dg1_lim_td[0])/dg1_lim_td[0]
mmr_lim_td_diff = np.abs(mmr_lim_td - mmr_lim_td[0])/mmr_lim_td[0]

times = np.asarray(no_lim_nc['time'])

plt.figure()
plt.semilogy(times, no_lim_td_diff, label='No limiter')
plt.semilogy(times, dg1_lim_td_diff, label='Standard limiter')
plt.semilogy(times, mmr_lim_td_diff, label='MMR limiter')
plt.ylabel('Difference in tracer density')
plt.xlabel('Time (s)')
plt.legend()

savename = f'{figure_stem}mmr_test_tracer_density5.jpg'
plt.savefig(savename, bbox_inches='tight')
print(f'saved figure to {savename}')

plt.close()