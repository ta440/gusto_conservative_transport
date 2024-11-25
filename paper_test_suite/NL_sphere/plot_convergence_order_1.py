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

order = 1
configuration = 'convergence'
quantity = 'm_X' # tracer_density or m_X

# Tracer density will be the total and m_X is the L2 of stead state error

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
if quantity == 'm_X':
    plot_name = f'{plot_dir}/NL_slice_ord_1_conv_m_l2_sse.png'
else:
    plot_name = f'{plot_dir}/NL_slice_ord_1_conv_Td_diff.png'
conservative = []
advective = []

# Spatial refinement levels (cells_per_edge = 2^ref_lev)
ncells_per_edge = [8, 16, 32]
refinement = [1/val for val in ncells_per_edge]

# Make a line for the conservative transport of the mixing ratio
print('Extracting conservative data')
for cell_no in ncells_per_edge:
    dirname = 'NL_sphere_conservative_order_'+str(order)+'_'+configuration+'_ncells_'+str(cell_no)
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
        print(np.log10(np.abs(m_err_l2[-1])))
    else:
        raise ValueError('Not correct inputs')

# Make a line for the advective transport of the mixing ratio
print('Extracting advective data')
for cell_no in ncells_per_edge:
    dirname = 'NL_sphere_advective_order_'+str(order)+'_'+configuration+'_ncells_'+str(cell_no)
    
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
        print(np.log10(np.abs(m_err_l2[-1])))
    else:
        raise ValueError('Not correct inputs')
    
##################################
# Labels for the two plots:
all_error_data = [advective, conservative]
colours = ['red', 'blue']
markers = ['s', 'o']
labels = ['advective', 'tracer conservative']

print(all_error_data)

xlabel = r"$log(\Delta x)$"
if quantity == 'tracer_density':
    log_by = 'axes'
    log_base=10
    ylabel = r"$log((T_d(T_{end}) - T_d(0))/T_d(0))$"
else:
    log_by = 'data'
    log_base= 'e'
    ylabel = "$log(||m(T_{end}) - m(0)||/||m(0)||$"


set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, refinement, error_data, label=label,
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
