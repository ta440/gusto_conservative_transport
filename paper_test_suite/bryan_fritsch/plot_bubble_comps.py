"""
Plots the moist Bryan & Fritsch bubble test case.

Here, we compare the conservatively and 
non-conservatively transported bubbles side-by-side.
"""
from os.path import abspath, dirname
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_fig, tomplot_field_title, extract_gusto_coords,
    extract_gusto_field
)

test1 = 'bryan_fritsch_advective_order_1'
test2 = 'bryan_fritsch_conservative_order_1'

results_file_name1 = f'{abspath(dirname(__file__))}/results/{test1}/field_output.nc'
results_file_name2 = f'{abspath(dirname(__file__))}/results/{test2}/field_output.nc'

plot_dir = f'{abspath(dirname(__file__))}/figures'

# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #
# Time index to compare at
t_idx = 2.

field_names = ['Theta_e', 'Theta_e']
time_idxs = [t_idx, t_idx]
cbars = [True, True]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contours = np.linspace(315.0, 325.0, 21)
remove_contour = 320.0
colour_scheme = 'RdBu_r'
field_label = r'$\theta_e$ (K)'
contour_method = 'tricontour'
xlims = [0., 10.]
ylims = [0., 10.]

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #

data_file1 = Dataset(results_file_name1, 'r')
data_file2 = Dataset(results_file_name2, 'r')

fig, axarray = plt.subplots(1, 2, figsize=(15, 6), sharex='all', sharey='all')

for i, (ax, time_idx, field_name, cbar) in \
        enumerate(zip(axarray.flatten(), time_idxs, field_names, cbars)):
    print(i)
    
    if i == 0:
        data = data_file1
    elif i == 1:
        data = data_file2

    # Data extraction ----------------------------------------------------------
    field_data = extract_gusto_field(data, field_name, time_idx=time_idx)
    coords_X, coords_Y = extract_gusto_coords(data, field_name)
    time = data['time'][time_idx]

    # Plot data ----------------------------------------------------------------
    cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=remove_contour)
    cf, lines = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines
    )

    if cbar:
        add_colorbar_fig(
            fig, cf, field_label, ax_idxs=[i], location='right'
        )
    tomplot_field_title(
        ax, f't = {time:.1f} s', minmax=True, field_data=field_data
    )

    # Labels -------------------------------------------------------------------
    if i == 0:
        ax.set_ylabel(r'$z$ (km)', labelpad=-20)
        ax.set_ylim(ylims)
        ax.set_yticks(ylims)
        ax.set_yticklabels(ylims)

    ax.set_xlabel(r'$x$ (km)', labelpad=-10)
    ax.set_xlim(xlims)
    ax.set_xticks(xlims)
    ax.set_xticklabels(xlims)

# Save figure ------------------------------------------------------------------

plot_name= f'{plot_dir}/bubble_comp_order'+str(order)+'_t'+str(time)+'s.png'
fig.subplots_adjust(wspace=0.15)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()