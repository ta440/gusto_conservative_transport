'''
Compare fields from the MMR test.
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


no_lim_results_dir = 'mmr_test_no_limiter_vary_rho'
dg1_lim_results_dir = 'mmr_test_standard_limiter_vary_rho'
mmr_results_dir = 'mmr_test_mmr_limiter_vary_rho'

ext_name = ''

no_lim_file_name = f'{abspath(dirname(__file__))}/results/{no_lim_results_dir}/field_output.nc'
dg1_lim_file_name = f'{abspath(dirname(__file__))}/results/{dg1_lim_results_dir}/field_output.nc'
mmr_lim_file_name = f'{abspath(dirname(__file__))}/results/{mmr_results_dir}/field_output.nc'

figure_stem = f'{abspath(dirname(__file__))}/figures/'

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
no_lim_file = Dataset(no_lim_file_name, 'r')
dg1_lim_file = Dataset(dg1_lim_file_name, 'r')
mmr_lim_file = Dataset(mmr_lim_file_name, 'r')

t_idxs = [0,-1,-1,-1]
title_names = ['Initial field', 'No limiter', 'DG1 limiter', 'MMR limiter']

colour_scheme = 'PiYG'
contour_method = 'contour'  # Need to use this method to show mountains!
field_name = 'm_X'

fig, axarray = plt.subplots(2,2, figsize=(8,8), constrained_layout='True', sharex='all', sharey='all')

for i, (ax, t_idx, title_name) in enumerate(zip(axarray.flatten(), t_idxs, title_names)):
    
    
    if i < 2:
        data = extract_gusto_field(no_lim_file, field_name, time_idx=t_idx)
    elif i == 2:
        data = extract_gusto_field(dg1_lim_file, field_name, time_idx=t_idx)
    elif i == 3:
        data = extract_gusto_field(mmr_lim_file, field_name, time_idx=t_idx)

    coords_X, coords_Y = extract_gusto_coords(no_lim_file, field_name, units = 'm')
    field_data, coords_X, coords_Y = \
        reshape_gusto_data(data, coords_X, coords_Y)


    #contours = tomplot_contours(field_data)
    contours = np.linspace(0.0,0.1,11)

    print(contours)

    cmap, lines = tomplot_cmap(contours, colour_scheme, extend_cmap='min')

    # Plot data ----------------------------------------------------------------
    cf, _ = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, 
    )   
    if i == 3:
        add_colorbar_ax(ax, cf, None, location='right')
    tomplot_field_title(ax, f'{title_name} \n', minmax=True, minmax_format='.4f', field_data=field_data)

    ax.set_aspect('equal')

    if i in [0, 2]:
        ax.set_ylabel(r'$z$', labelpad=0)
        #ax.set_yticks([0,0.025,0.05,0.075,0.1])
        #ax.set_yticklabels([0,0.025,0.05,0.075,0.1])


savename = f'{figure_stem}mmr_test_compare_{field_name}.jpg'
plt.savefig(savename, bbox_inches='tight')
print(f'saved figure to {savename}')

plt.close()