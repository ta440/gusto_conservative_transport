# A script to make a convergence plot from 
# the conservation test with the NL_slice test.
# This uses the tomplot library

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

order = 0
configuration = 'convergence'
quantity = 'm_X' # tracer_density or m_X

# Tracer density will be the total and m_X is the L2 of stead state error

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
if quantity == 'm_X':
    if configuration == 'convergence':
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_conv_m_l2sse.png'
    else:
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_consist_m_l2sse.png'
else:
    if configuration == 'convergence':
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_conv_Td_diff.png'
    else:
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_consist_Td_diff.png'

# Arrays to store plotting values
conservative = []
advective = []

# Spatial refinements
dxzs = [120, 140, 160, 180, 200]# [100, 140, 160, 180]
dx_values = []

for dxz in dxzs:
    dx_values.append(2000./dxz)

# Make a line for the conservative transport of the mixing ratio
for dxz in dxzs:
    dirname = 'NL_slice_conservative_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    nc = Dataset(f'results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]
        
        # Compute the relative change in rho_X:
        T_d_diff = np.abs(T_d[-1] - T_d[0])/T_d[0]
        
        conservative.append(T_d_diff)
    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        conservative.append(np.abs(m_err_l2[-1]))
    else:
        raise ValueError('Not correct inputs')

# Make a line for the advective transport of the mixing ratio
for dxz in dxzs:
    dirname = 'NL_slice_advective_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    
    nc = Dataset(f'results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]
        
        # Compute the relative change in rho_X:
        T_d_diff = np.abs(T_d[-1] - T_d[0])/T_d[0]
        
        advective.append(T_d_diff)
    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        advective.append(np.abs(m_err_l2[-1]))
    else:
        raise ValueError('Not correct inputs')
    
##################################
# Labels for the two plots:
all_error_data = [advective, conservative]
colours = ['red', 'blue']
markers = ['s', 'o']
labels = ['advective', 'tracer conservative']

print(all_error_data)

if quantity == 'tracer_density':
    log_by = 'axes'
    log_base=10
    ylabel = r"$(T_d(T_{end}) - T_d(0)/T_d(0))$"
    xlabel = r"$log(\Delta x)$"
else:
    log_by='data'
    log_base='e'
    ylabel = "$ln(||m(T_{end}) - m(0)||/||m(0)||$"
    xlabel = r"$ln(\Delta x)$"

set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, dx_values, error_data, label=label,
                     color=colour, marker=marker, log_by=log_by,
                     log_base=log_base)
                     
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
