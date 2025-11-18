"""
Plot from the Terminator Toy test case.
We examine the distribution of the dry density (rho_d)
 and species (X and X2) at different times.
This plots:
(a) rho_d @ t = 0 days, (b) rho_d @ t = 518400 s (6 days), (c) rho_d @ t = 1036800 s (12 days, final time)
(d) X @ t = 0 days, (e) X @ t = 518400 s (6 days), (f) X @ t = 1036800 s (12 days, final time)
(g) X2 @ t = 0 days, (h) X2 @ t = 518400 s (6 days), (i) X2 @ t = 1036800 s (12 days, final time)
"""
from os.path import abspath, dirname
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field, add_colorbar_fig,
    tomplot_field_title, extract_gusto_coords, extract_gusto_field
)
import cartopy.crs as ccrs

tracer_spaces = 'same'
prog = 'rho_d'
ncells_per_edge = 24


# Not with consistent projection
#adv_extra_name = 'analyt_ref_Td_solve_'
#con_extra_name = 'aug1_zerocrops_and_chemlim'
#adv_file_name = f'/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/terminator_toy/results/terminator_toy_{adv_extra_name}advective_ncells_{ncells_per_edge}/field_output.nc'
#con_file_name = f'/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/terminator_toy/results/terminator_toy_{con_extra_name}conservative_ncells_{ncells_per_edge}/field_output.nc'
#extra_name = ''

# Now with consistent projection
adv_extra_name = 'analyt_ref_Td_solve_'
con_extra_name = 'nov5_consistent_proj_'
adv_file_name = f'/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/terminator_toy/results/terminator_toy_{adv_extra_name}advective_ncells_{ncells_per_edge}/field_output.nc'
con_file_name = f'{abspath(dirname(__file__))}/results/terminator_toy_{con_extra_name}conservative_ncells_{ncells_per_edge}/field_output.nc'
extra_name = '_nov5_consistent_proj'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these paths need editing, which will usually involve
# removing the abspath part to set directory paths relative to this file
#results_file_name = f'{abspath(dirname(__file__))}/results/terminator_toy_{tracer_spaces}_spaces_{scheme}_ncells_{ncells_per_edge}/field_output.nc'
#plot_stem = f'{abspath(dirname(__file__))}/figures/terminator_toy_all_{tracer_spaces}_spaces_{scheme}_ncells_{ncells_per_edge}'

plot_stem = f'{abspath(dirname(__file__))}/figures/terminator_toy_robinson_comparison{extra_name}'


# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #

field_names = ['rho_d', 'rho_d', 'rho_d',
               'X_tracer', 'X_tracer', 'X_tracer',
               'X2_tracer', 'X2_tracer', 'X2_tracer']
time_idxs = [0, -1, -1,
             0, -1, -1,
             0, -1, -1]
cbars = [False, False, True,
         False, False, True,
         False, False, True]
schemes= ['adv','adv','con',
         'adv','adv','con',
         'adv','adv','con']

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #

projection = ccrs.Robinson()

# Save filled contour data for adding colourbars at end
all_cf = []

rho_contours = np.linspace(0.0, 1.50, 15)
rho_colour_scheme = 'gist_earth_r'
rho_field_label = r'$\rho$ (kg m$^{-3}$)'

X_contours_first = np.linspace(0.0, 4.009e-6, 15)
X_contours = np.linspace(0.0, 4.0e-6, 15)
X_colour_scheme = 'Purples'
X_field_label = r'$X$ (kg kg$^{-1}$)'

X2_contours_first = np.linspace(0.0, 2.000000000001e-6, 15)
X2_contours = np.linspace(0.0, 2.0e-6, 15)
X2_colour_scheme = 'Greens'
X2_field_label = r'$X_2$ (kg kg$^{-1}$)'

contour_method = 'tricontour'
xlims = [-180, 180]
ylims = [-90, 90]

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
adv_data_file = Dataset(adv_file_name, 'r')
con_data_file = Dataset(con_file_name, 'r')

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #
subplots_x = 3
subplots_y = 3
fig = plt.figure(figsize=(20, 12))

for i, (time_idx, field_name, cbar, scheme) in \
        enumerate(zip(time_idxs, field_names, cbars, schemes)):

    ax = fig.add_subplot(subplots_y, subplots_x, 1+i, projection=projection)
    
    if scheme == 'adv':
        data_file = adv_data_file
    else:
        data_file = con_data_file

    if time_idx == 'midpoint':
        time_idx = int((len(data_file['time'][:]) - 1) / 2)

    # Data extraction ----------------------------------------------------------
    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)

    print(min(field_data))
    coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
    # Quote time in days:
    time = data_file['time'][time_idx] / (24.*60.*60.)

    # Select options for each field --------------------------------------------
    if field_name == 'rho_d':
        contours = rho_contours
        colour_scheme = rho_colour_scheme
        field_label = rho_field_label
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None, extend_cmap=True)
        cbar_labelpad = -80
        data_format = '.2e'

    elif field_name == 'X_tracer':
        if time_idx == 0:
            contours = X_contours_first
        else:
            contours = X_contours
        colour_scheme = X_colour_scheme
        field_label = X_field_label
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None, extend_cmap=True)
        cbar_labelpad = -80
        data_format = '.2e'

    elif field_name == 'X2_tracer':
        if time_idx == 0:
            contours = X2_contours_first
        else:
            contours = X2_contours
        colour_scheme = X2_colour_scheme
        field_label = X2_field_label
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None, extend_cmap=True)
        cbar_labelpad = -80
        data_format = '.2e'

    # Plot data ----------------------------------------------------------------

    cf, _ = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, projection=projection
    )

    all_cf.append(cf)

    if i == 0:
        title = 'Initial condition \n '
    elif i == 1:
        title = f'Advective scheme, t = {time:.1f} days \n '
    elif i == 2:
        title = f'Tracer conservative scheme, t = {time:.1f} days \n '
    else:
        title = ' '
    tomplot_field_title(
        ax, title, minmax=True, minmax_format=data_format,
        field_data=field_data
    )

    # Labels -------------------------------------------------------------------
    #if i in [0, 3, 6]:
    #    ax.set_ylabel(r'$\vartheta$ (deg)', labelpad=-20)
    #    ax.set_ylim(ylims)
    #    ax.set_yticks(ylims)
    #    ax.set_yticklabels(ylims)

    #if i in [6, 7, 8]:
    #    ax.set_xlabel(r'$\lambda$ (deg)', labelpad=-10)
    #    ax.set_xlim(xlims)
    #    ax.set_xticks(xlims)
    #    ax.set_xticklabels(xlims)

for i, (cbar, field_name) in enumerate(zip(cbars, field_names)):

    cbar_labelpad = -50
    data_format = '.1e'

    # Get information for field
    if field_name == 'rho_d':
        field_label = rho_field_label

    elif field_name == 'X_tracer':
        field_label = X_field_label

    elif field_name == 'X2_tracer':
        field_label = X2_field_label

    if cbar:
        add_colorbar_fig(
            fig, all_cf[i], field_label, ax_idxs=[i-1, i], location='right',
            cbar_labelpad=cbar_labelpad, cbar_format=data_format
        )


# Save figure ------------------------------------------------------------------
fig.subplots_adjust(wspace=0.25)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()