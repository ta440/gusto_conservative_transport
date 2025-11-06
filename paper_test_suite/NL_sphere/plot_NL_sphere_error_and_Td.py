# A script to make a convergence plot from 
# the conservation test with the NL_slice test.
# This uses the tomplot library
# This is for the convergence configuration of the tests.

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig, add_convergence_comparison_line)
import numpy as np

order = 1
quantity = 'm_X' # tracer_density or m_X

# Tracer density will be the total and m_X is the L2 of stead state error

# Older branch:
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/NL_sphere/'
results_dir = '/data/home/ta440/firedrake_pip_19062025/src/gusto/gusto_conservative_transport/paper_test_suite/NL_sphere/'
#name_ext=''
# Or, for current directory:
#results_dir = ''

# For any additional components to the name:
#name_ext = '_dt_450.0_divflow_Vu_CG1'

# order = 0 local branch
#name_ext = '_dt_450.0_divflow_Vu_CG1_td_solve'

# extra names
if order == 0:
    extra_name = 'proj_18aug_'
    name_ext = '_dt_450.0'
elif order == 1:
    extra_name = ''
    name_ext = '_td_solve'


# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
if quantity == 'm_X':
    plot_name = f'{plot_dir}/NL_sphere_ord_{order}_convergence_m_l2_sse.png'
else:
    plot_name = f'{plot_dir}/NL_sphere_ord_{order}_convergence_Td_sum.png'
conservative = []
advective = []

# Spatial resolutions.
# Order 1 can be twice as coarse as order 0.
if order == 1:
    ncells_per_edge = [8,16,24,32]
    #refinement = ncells_per_edge
    refinement = [32/val for val in ncells_per_edge]
else:
    #ncells_per_edge = [36,48]
    ncells_per_edge = [16,24,36,48]
    refinement = [48/val for val in ncells_per_edge]

# Make a line for the conservative transport of the mixing ratio
print('Extracting conservative data')
for cell_no in ncells_per_edge:
    dirname = f'NL_sphere_{extra_name}conservative_order_{order}_convergence_ncells_{cell_no}{name_ext}'
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

    if quantity == 'tracer_density':
    # Extract the tracer density:
    
        rho_X = nc.groups['TracerDensity_m_X_rho_d']
        
        T_d = rho_X.variables['total'][:]
        
        # Just compute Tracer density at the end time:
        #T_d_diff = np.abs(T_d[-1] - T_d[0])/T_d[0]
        #conservative.append(T_d_diff)

        # Compute Tracer Density as a sum:
        T_d_sum = np.sum(np.abs(T_d[1:] - T_d[0])/T_d[0])/len(T_d[1:])
        conservative.append(T_d_sum)

    elif quantity == 'm_X':
        m_err = nc.groups['m_X_error']
        m_err_l2 = m_err.variables['l2'][:]
        # Normalise steady state error by initial value
        # of m_X
        mX = nc.groups['m_X']
        mX0 = mX.variables['l2'][0]
        print(mX0)
        conservative.append(np.abs(m_err_l2[-1])/mX0)
        print(np.log(np.abs(m_err_l2[-1])))
    else:
        raise ValueError('Not correct inputs')

# Make a line for the advective transport of the mixing ratio
print('Extracting advective data')
for cell_no in ncells_per_edge:
    dirname = f'NL_sphere_{extra_name}advective_order_{order}_convergence_ncells_{cell_no}{name_ext}'
    
    nc = Dataset(f'{results_dir}results/{dirname}/diagnostics.nc')

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
    if order == 1:
        xlabel = r"$32/N_e$"
    else:
        xlabel = r"$48/N_e$"
    gradient_in_label = False
    labels = ['advective', 'tracer conservative']
else:
    log_by = 'data'
    log_base= 'e'
    if order == 1:
        xlabel = r"ln(32/$_e$)"
    else:
        xlabel = r"ln(48/$N_e$)"
    ylabel = r"ln(Final Mixing Ratio L2 error)"
    gradient_in_label = True
    labels = ['advective: ', 'tracer conservative: ']


set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, refinement, error_data, label=label,
                    color=colour, marker=marker, log_by=log_by,
                    log_base=log_base, gradient_in_label=gradient_in_label)
    
if quantity == 'm_X':   
    if order == 1:
        x_points = [1,4]
        y_shift = 0.5
    elif order == 0:
        x_points = [1,3]
        y_shift = -0.5
    add_convergence_comparison_line(ax, 2, label=r'$(\Delta x)^2$', color='k', log_by=log_by,
                                    x_points=x_points, y_shift=y_shift)
                    
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

#xtick_labels = [np.round(val,2) for val in refinement]
#xticks = [np.round(val,2) for val in refinement]
#ax.set_xticks(xticks, xtick_labels)

tomplot_legend_ax(ax, location='bottom')

plt.grid()

if quantity == 'tracer_density':
    old_xticklabels = ax.get_xticklabels()

    # Set most tick labels to be an empty string
    new_xticklabels = ['' for _ in old_xticklabels]

    ax.set_xticklabels([])
    #ax.set_xticklabels(new_xticklabels)

    xtick_labels = [np.round(val,2) for val in refinement]
    #print(xtick_labels)
    xticks = [np.round(val,2) for val in refinement]
    ax.set_xticks(xticks, xtick_labels)

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()
