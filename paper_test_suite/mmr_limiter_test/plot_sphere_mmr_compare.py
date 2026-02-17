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
    extract_gusto_coords, extract_gusto_field, reshape_gusto_data,
    add_colorbar_fig
)


no_lim_results_dir = 'mmr_test_no_limiter_sphere_cylinder'
mmr_results_dir = 'mmr_test_mmr_limiter_sphere_cylinder'

ext_name = ''

no_lim_file_name = f'{abspath(dirname(__file__))}/results/{no_lim_results_dir}/field_output.nc'
mmr_lim_file_name = f'{abspath(dirname(__file__))}/results/{mmr_results_dir}/field_output.nc'

figure_stem = f'{abspath(dirname(__file__))}/figures/'

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
no_lim_file = Dataset(no_lim_file_name, 'r')
mmr_lim_file = Dataset(mmr_lim_file_name, 'r')

t_idxs = [0,-1,-1]
title_names = ['Initial field', 'No limiter', 'MMR limiter']

colour_scheme = 'PiYG'
contour_method = 'tricontour'  

domain_limit = {'X' : (-180, 180), 'Y' : (-90, 90)}
y_ticks = [-90, 0, 90]
y_tick_labels = ['90S', '0', '90N']
xticks = [-180, 0, 180]
xtick_labels = ['180W', '0', '180E']

field_name = 'm_X'


fig, axarray = plt.subplots(1,3, figsize=(10,3), constrained_layout='True', sharey='all')

for i, (ax, t_idx, title_name) in enumerate(zip(axarray.flatten(), t_idxs, title_names)):
    
    
    if i < 2:
        data = extract_gusto_field(no_lim_file, field_name, time_idx=t_idx)
    elif i == 2:
        data = extract_gusto_field(mmr_lim_file, field_name, time_idx=t_idx)

    coords_X, coords_Y = extract_gusto_coords(no_lim_file, field_name)
    field_data, coords_X, coords_Y = reshape_gusto_data(data, coords_X, coords_Y)


    #contours = tomplot_contours(field_data)
    contours = np.linspace(0.0,1.2,7)
    cbar_ticks = contours

    print(contours)

    cmap, lines = tomplot_cmap(contours, colour_scheme, extend_cmap='min')

    # Plot data ----------------------------------------------------------------
    cf, _ = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, 
        )   
    #if i == 1:
    #    add_colorbar_ax(axarray, cf, r"$m$ (kg kg$^{-1}$)", location='bottom', cbar_ticks=cbar_ticks, cbar_format='.1f', pad=0.07)
    
    
    
    
    tomplot_field_title(ax, f'{title_name} \n', minmax=True, minmax_format='.3f', field_data=field_data)

    ax.set_aspect('equal')

    #if i == 0:
     #   ax.set_ylabel(r'lat (deg)', labelpad=0)
    #ax.set_xlabel(r'lon (deg)', labelpad=0)

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize = '15.0')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize='15.0')

add_colorbar_fig(fig, cf, r"$m$ (kg kg$^{-1}$)", location='bottom', cbar_ticks=cbar_ticks, cbar_format='.1f', 
                 cbar_padding=0.15, cbar_thickness=0.05)

savename = f'{figure_stem}mmr_sphere_compare_{field_name}.jpg'
plt.savefig(savename, bbox_inches='tight')
print(f'saved figure to {savename}')

plt.close()