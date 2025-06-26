# A script to make a convergence plot from 
# the conservation test with the NL_slice test.
# This uses the tomplot library
# This is to plot the consistency configuration over time.

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

order = 1


# Tracer density will be the total and m_X is the L2 of stead state error

# Older branch:
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/NL_sphere/'

# Or, for current directory:
results_dir = ''

# For any additional components to the name:
extra_name = ''

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_name = f'{abspath(dirname(__file__))}/figures/NL_sphere_ord_{order}_consistency_m_l2_sse_over_time.png'

if order == 0:
    Ne = 48
elif order == 1:
    Ne = 32

# Advective result:
adv_dirname = f'NL_sphere_advective_order_{order}_consistency_ncells_{Ne}{extra_name}'
nc = Dataset(f'{abspath(dirname(__file__))}/results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])
time_days = times/60/60/24

adv_mx_sse = nc.groups['m_X_error']
adv_mx_sse = adv_mx_sse['l2'][:]

# Normalise:
adv_mX = nc.groups['m_X']
adv_mX = adv_mX['l2'][:]
adv_mX0 = adv_mX[0]
adv_mx_sse = adv_mx_sse/adv_mX0


# Conservative result:
con_dirname = f'NL_sphere_conservative_order_{order}_consistency_ncells_{Ne}{extra_name}'
conservative_nc = Dataset(f'{abspath(dirname(__file__))}/results/{con_dirname}/diagnostics.nc')

con_mx_sse = conservative_nc.groups['m_X_error']
con_mx_sse = con_mx_sse['l2'][:]

# Normalise:
con_mX = nc.groups['m_X']
con_mX = con_mX['l2'][:]
con_mX0 = con_mX[0]
con_mx_sse = con_mx_sse/con_mX0


# Plot these tracer densities
plt.figure()
plt.semilogy(time_days, adv_mx_sse, label='advective', c='r')
plt.semilogy(time_days, con_mx_sse, label='tracer conservative', c='b')
plt.xlabel('Time (days)', size=12)
plt.ylabel('Relative mixing ratio error', size=12)
plt.legend(loc='lower right', prop={'size': 12})

print(f'Saving figure to {plot_name}')
plt.savefig(plot_name, bbox_inches='tight')
plt.close()