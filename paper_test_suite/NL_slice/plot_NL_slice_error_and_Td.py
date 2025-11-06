# A script to make a convergence plot from 
# the conservation test with the NL_slice test.
# This uses the tomplot library
# This script is currently being used to generate the best 
# plots for the paper.

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig, add_convergence_comparison_line)
import numpy as np
from matplotlib.ticker import ScalarFormatter, NullLocator

order = 0
configuration = 'convergence'
quantity = 'm_X' # tracer_density or m_X

#Either link to current results or from previous branch
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/NL_slice'
results_dir = '/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/NL_slice'

#extra_name = 'Tdsolve'
extra_name_adv = 'proj_18aug_'
extra_name_con = 'proj_18aug_'

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
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_conv_Td_sum.png'
    else:
        plot_name = f'{plot_dir}/NL_slice_ord_{order}_consist_Td_sum.png'

# Arrays to store plotting values
conservative = []
advective = []

# Spatial refinements
if order == 1:
    dxzs = [50, 60, 70, 80, 90, 100]
else:
    dxzs = [100, 120, 140, 160, 180, 200]

dx_values = []

for dxz in dxzs:
    dx_values.append(2000./dxz)

# Make a line for the conservative transport of the mixing ratio
for dxz in dxzs:
    dirname = f'NL_slice_{extra_name_con}conservative_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    nc = Dataset(f'{results_dir}/results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]
        
       # Just compute Tracer density at the end:
        #T_d_diff = np.abs(T_d[-1] - T_d[0])/T_d[0]
        #conservative.append(T_d_diff)

        # Compute Tracer Density as a sum:
        T_d_sum = np.sum(np.abs(T_d[1:] - T_d[0])/T_d[0])/len(T_d[1:])
        conservative.append(T_d_sum)

    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        mX = nc.groups['m_X']
        mX0 = mX.variables['l2'][0]
        conservative.append(np.abs(m_err_l2[-1])/mX0)
    else:
        raise ValueError('Not correct inputs')

# Make a line for the advective transport of the mixing ratio
for dxz in dxzs:
    dirname = f'NL_slice_{extra_name_adv}advective_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    
    nc = Dataset(f'{results_dir}/results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]
        
        # Just compute Tracer density at the end:
        #T_d_diff = np.abs(T_d[-1] - T_d[0])/T_d[0]
        #advective.append(T_d_diff)

        # Compute Tracer Density as a sum:
        T_d_sum = np.sum(np.abs(T_d[1:] - T_d[0])/T_d[0])/len(T_d[1:])
        advective.append(T_d_sum)
    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        mX = nc.groups['m_X']
        mX0 = mX.variables['l2'][0]
        print(mX0)
        advective.append(np.abs(m_err_l2[-1])/mX0)
        print(np.log(np.abs(m_err_l2[-1])))
    else:
        raise ValueError('Not correct inputs')
    
##################################
# Labels for the two plots:
all_error_data = [advective, conservative]
colours = ['red', 'blue']
markers = ['s', 'o']

print(all_error_data)

if quantity == 'tracer_density':
    log_by = 'axes'
    log_base = 10
    #ylabel = r"$T_d(T_{end}) - T_d(0)/T_d(0)$"
    ylabel = 'Mean Tracer Density Error'
    xlabel = r"$\Delta x$ (m)"
    gradient_in_label = False
    labels = ['advective', 'tracer conservative']
else:
    log_by='data'
    log_base='e'
    ylabel = r"ln(Final Mixing Ratio L2 error)"
    xlabel = r"ln$(\Delta x)$ (ln(m))"
    gradient_in_label = True
    labels = ['advective: ', 'tracer conservative: ']

set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, dx_values, error_data, label=label,
                     color=colour, marker=marker, log_by=log_by,
                     log_base=log_base, gradient_in_label=gradient_in_label)
                     
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

if quantity == 'm_X':   
    if order == 1:
        x_points = [20,40]
        y_shift = 0.5
    elif order == 0:
        x_points = [10,20]
        y_shift = 0.2
    add_convergence_comparison_line(ax, 2, label=r'$(\Delta x)^2$', color='k', log_by=log_by,
                                    x_points=x_points, y_shift=y_shift)

if order == 1:  
    if quantity == 'tracer_density':
        ax.set_xlim([19,41])
        ax.set_ylim([1e-14,1e-4])
        plt.gca().axes.xaxis.set_ticks([20,25,30,35,40])
        plt.gca().axes.xaxis.set_ticklabels([20,25,30,35,40])
        plt.gca().xaxis.set_minor_formatter(ScalarFormatter())
else:  
    if quantity == 'tracer_density':
        ax.set_xlim([9.5, 20.5])
        ax.set_ylim([1e-13,1e-4])
        #plt.gca().xaxis.set_minor_formatter(ScalarFormatter())
        #plt.gca().xaxis.set_major_formatter(ScalarFormatter())
        #plt.xticks([],minor=False)
        #plt.xticks([10,12,14,16], minor=False)
        plt.gca().axes.xaxis.set_ticks([10,12,14,16,18,20])
        plt.gca().axes.xaxis.set_ticklabels([10,12,14,16,18,20])
        plt.gca().xaxis.set_minor_locator(NullLocator())

tomplot_legend_ax(ax, location='bottom')

plt.grid()
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()
