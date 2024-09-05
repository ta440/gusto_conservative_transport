# A script to make a convergence plot from 
# the conservation test with the NL_slice test.
# This uses the tomplot library

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)

order = 1
configuration = 'convergence'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/NL_slice_order_1_convergence.png'

conservative = []
advective = []

# Spatial refinements
dxzs = [100, 140]
dx_values = [2000./100, 2000./140]

# Make a line for the conservative transport of the mixing ratio
for dxz in dxzs:
    dirname = 'NL_slice_conservative_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    nc = Dataset(f'results/{dirname}/diagnostics.nc')

    # Extract the tracer density:
    rho_X = nc.groups['TracerDensity_m_X_rho_d']
    
    T_d = rho_X.variables['total'][:]
    
    # Compute the relative change in rho_X:
    T_d_diff = (T_d[-1] - T_d[0])/T_d[0]
    
    conservative.append(T_d_diff)


# Make a line for the advective transport of the mixing ratio
for dxz in dxzs:
    dirname = 'NL_slice_advective_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    
    nc = Dataset(f'results/{dirname}/diagnostics.nc')

    # Extract the tracer density:
    rho_X = nc.groups['TracerDensity_m_X_rho_d']

    T_d = rho_X.variables['total'][:]
    
    # Compute the relative change in rho_X:
    T_d_diff = (T_d[-1] - T_d[0])/T_d[0]
    
    advective.append(T_d_diff)
    
##################################
# Labels for the two plots:
all_error_data = [advective, conservative]
colours = ['red', 'blue']
markers = ['s', 'o']
labels = ['advective', 'tracer conservative']

print(all_error_data)

log_by = 'data'
xlabel = r"$\log(\Delta x)$"
ylabel = r"$\log(|T_d(T_{end}) - T_d(0)|/|T_d(0)|)$"
    
set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, dx_values, error_data, label=label,
                     color=colour, marker=marker, log_by=log_by)
                     
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

tomplot_legend_ax(ax, location='bottom')

plt.grid()
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()
