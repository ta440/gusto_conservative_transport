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

order = 0

test1 = f'bryan_fritsch_linear_mX0_16x1_Tdsolve_advective_order_{order}dxz50'
test2 = f'bryan_fritsch_linear_mX0_16x1_Tdsolve_conservative_order_{order}dxz50'

results_dir = '/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/bryan_fritsch/'

results_file_name1 = f'{results_dir}results/{test1}/field_output.nc'
results_file_name2 = f'{results_dir}results/{test2}/field_output.nc'

plot_dir = f'{abspath(dirname(__file__))}/figures'

# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #
# Time index to compare at
t_idx = 10.
field = 'Theta_e'

field_names = [field, field]
time_idxs = [t_idx, t_idx]
cbars = [False, True]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
if field == 'theta':
    contours = np.linspace(290,330,21)
elif field == 'rho':
    contours = np.linspace(0.4,1.2,21)
elif field == 'Theta_e':
    contours = np.linspace(317.0, 326.0, 9)
elif field == 'cloud_water':
    contours = np.linspace(0.01,0.02,21)
elif field == 'water_vapour':
    contours = np.linspace(0.0,0.01,21)
remove_contour = None#320.0
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
        cbar_vals = np.linspace(contours[0], contours[-1], 5)
        add_colorbar_fig(
            fig, cf, field_label, ax_idxs=[i], location='right',
            cbar_ticks=cbar_vals, cbar_format='.0f'
        )

    if i ==0:
        tomplot_field_title(
            ax, f'Advective transport \n',  minmax=True, field_data=field_data
        )
    else:
        tomplot_field_title(
            ax, f'Conservative transport \n', minmax=True, field_data=field_data
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

plot_name= f'{plot_dir}/bubble_comp_order'+str(order)+'_t'+str(time)+'s_'+field+'.png'
fig.subplots_adjust(wspace=0.15)
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()