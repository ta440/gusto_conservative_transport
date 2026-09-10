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
#quantity = 'tracer_density'
quantity = 'm_X'

# Tracer density will be the total and m_X is the L2 of stead state error

# Older branch:
adv_results_dir = '/data/home/ta440/firedrake_pip_29_1_26/gusto/gusto_conservative_transport/paper_test_suite/NL_sphere/'
con_results_dir = ''

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
    adv_name_ext = ''
    con_name_ext=''


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

# Spatial resolutions
if order == 1:
    ncells_per_edge = [16,24,32,48]
    #refinement = [32/val for val in ncells_per_edge]
else:
    #refinement = [48/val for val in ncells_per_edge]
    ncells_per_edge = [16,24,36,48]

R = 6371220. 
refinement = [np.pi*R/(2*Ne)/1000 for Ne in ncells_per_edge]
print(refinement)

# Make a line for the conservative transport of the mixing ratio
print('Extracting conservative data')
for cell_no in ncells_per_edge:
    dirname = f'NL_sphere_{extra_name}conservative_order_{order}_convergence_ncells_{cell_no}{con_name_ext}'
    nc = Dataset(f'{con_results_dir}results/{dirname}/diagnostics.nc')

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
    dirname = f'NL_sphere_{extra_name}advective_order_{order}_convergence_ncells_{cell_no}{adv_name_ext}'
    
    nc = Dataset(f'{adv_results_dir}results/{dirname}/diagnostics.nc')

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
    ylabel = r"$\overline{\Delta M}$"
    #ylabel = 'Mean Tracer Density Error'
    if order == 1:
        xlabel = r"$\Delta x$ (km)"
    else:
        xlabel = r"$\Delta x$ (km)"
    gradient_in_label = False
    labels = ['advective', 'conservative']
else:
    log_by = 'data'
    log_base= 'e'
    if order == 1:
        xlabel = r"ln($\Delta x$) (ln(km))"
    else:
        xlabel = r"ln($\Delta x$) (ln(km))"
    ylabel = r"ln($||m(T_{end}) - m(0)||_2/||m(0)||_2$)"
    #ylabel = r"ln(Final Mixing Ratio L2 error)"
    gradient_in_label = True
    labels = ['advective: ', 'conservative: ']



set_tomplot_style()
fig, ax = plt.subplots(1, 1, figsize=(5, 5)) 
for error_data, colour, marker, label in \
        zip(all_error_data, colours, markers, labels):
    plot_convergence(ax, refinement, error_data, label=label,
                    color=colour, marker=marker, log_by=log_by,
                    log_base=log_base, gradient_in_label=gradient_in_label)
    
if quantity == 'm_X':   
    if order == 1:
        x_points = None#[10,14]
        y_shift = 0.5
    elif order == 0:
        x_points = None#[1,3]
        y_shift = 0.5
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

    # Set most tick labels to be an empty string#
    new_xticklabels = ['' for _ in old_xticklabels]
    print(new_xticklabels)

    ax.set_xticks([])
    ax.set_xticklabels([])
    #ax.set_xticks([])

    if order == 1:
        ax.minorticks_off()
        #xtick_labels = [300, 600, 900, 1200]
        #xticks = [300, 600, 900, 1200]
        xtick_labels = [200, 300, 400, 500, 600, 700]
        xticks = [200, 300, 400, 500, 600, 700]
    else:
        xtick_labels = [200, 300, 400, 500, 600]
        xticks = [200, 300, 400, 500, 600]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    ax.set_ylim([1e-16,1e-3])

    #ax.xaxis.set_minor_formatter(plt.NullFormatter())


# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, dpi=500, bbox_inches='tight')
plt.close()
