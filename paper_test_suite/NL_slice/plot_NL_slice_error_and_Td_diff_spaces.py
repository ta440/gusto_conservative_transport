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

order = 1
configuration = 'convergence'
quantity = 'tracer_density' # tracer_density or m_X

# Results from older branch:
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/NL_slice/'

# Results in this branch:
results_dir = ''

extra_name_same = 'proj_18aug'
extra_name_diff = 'diff_spaces'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
if quantity == 'm_X':
    plot_name = f'{plot_dir}/NL_slice_{extra_name_diff}_ord_{order}_conv_m_l2sse.png'
else:
    plot_name = f'{plot_dir}/NL_slice_{extra_name_diff}_ord_{order}_conv_Td_sum.png'

# Arrays to store plotting values
conservative = []
advective = []
conservative_diff = []
advective_diff = []

# Spatial refinements
dxzs = [50, 60, 70, 80, 90, 100]

dx_values = []

for dxz in dxzs:
    dx_values.append(2000./dxz)

# Make a line for the conservative transport of the mixing ratio
for dxz in dxzs:
    #dirname = 'NL_slice_conservative_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    dirname = f'NL_slice_{extra_name_same}_conservative_order_{order}_{configuration}_dxz_{dxz}'
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]

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
    #dirname = 'NL_slice_advective_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    dirname = f'NL_slice_{extra_name_same}_advective_order_{order}_{configuration}_dxz_{dxz}'
    
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]

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
    

for dxz in dxzs:
    #dirname = 'NL_slice_conservative_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    dirname = f'NL_slice_{extra_name_diff}_conservative_order_{order}_{configuration}_dxz_{dxz}'
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]

        # Compute Tracer Density as a sum:
        T_d_sum = np.sum(np.abs(T_d[1:] - T_d[0])/T_d[0])/len(T_d[1:])
        conservative_diff.append(T_d_sum)

    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        mX = nc.groups['m_X']
        mX0 = mX.variables['l2'][0]
        conservative_diff.append(np.abs(m_err_l2[-1])/mX0)
    else:
        raise ValueError('Not correct inputs')

# Make a line for the advective transport of the mixing ratio
for dxz in dxzs:
    #dirname = 'NL_slice_advective_order_'+str(order)+'_'+configuration+'_dxz_'+str(dxz)
    dirname = f'NL_slice_{extra_name_diff}_advective_order_{order}_{configuration}_dxz_{dxz}'
    
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]

        # Compute Tracer Density as a sum:
        T_d_sum = np.sum(np.abs(T_d[1:] - T_d[0])/T_d[0])/len(T_d[1:])
        advective_diff.append(T_d_sum)
    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        mX = nc.groups['m_X']
        mX0 = mX.variables['l2'][0]
        print(mX0)
        advective_diff.append(np.abs(m_err_l2[-1])/mX0)
        print(np.log(np.abs(m_err_l2[-1])))
    else:
        raise ValueError('Not correct inputs')
    
##################################
# Labels for the two plots:
all_error_data = [advective, conservative, advective_diff, conservative_diff]
colours = ['red', 'blue', 'green', 'orange']
markers = ['s', 'o', 'x', 'v']
labels = ['advective (same space)', 'tracer conservative (same space)', \
          'advective (diff space)', 'tracer conservative (diff space)']

print(all_error_data)

if quantity == 'tracer_density':
    log_by = 'axes'
    log_base = 10
    #ylabel = r"$T_d(T_{end}) - T_d(0)/T_d(0)$"
    ylabel = 'Mean Tracer Density Error'
    xlabel = r"$\Delta x$ (m)"
else:
    log_by='data'
    log_base='e'
    #ylabel = "$ln(||m(T_{end}) - m(0)||/||m(0)||)$"
    ylabel = "$ln$(Final Mixing Ratio L2 error)"
    xlabel = r"$ln(\Delta x)$ ($ln$(m))"

set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, dx_values, error_data, label=label,
                     color=colour, marker=marker, log_by=log_by,
                     log_base=log_base, gradient_in_label=False)
                     
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

if order == 1:  
    if quantity == 'tracer_density':
        ax.set_xlim([19,41])
        ax.set_ylim([1e-15,1e-4])
        plt.gca().axes.xaxis.set_ticks([20,25,30,35,40])
        plt.gca().axes.xaxis.set_ticklabels([20,25,30,35,40])
        plt.gca().xaxis.set_minor_formatter(ScalarFormatter())
    else:
        add_convergence_comparison_line(ax, 2, label=r'$(\Delta x)^2$', color='black',
                                    log_by=log_by, log_base=log_base)
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
