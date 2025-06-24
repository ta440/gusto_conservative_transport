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
Ne = 32

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
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{abspath(dirname(__file__))}/NL_sphere_ord_{order}_consistency_m_l2_sse_over_time.png'

# Advective result:
adv_dirname = f'NL_sphere_advective_order_{order}_consistency_ncells_{Ne}{extra_name}/field_output.nc'
nc = Dataset(f'{abspath(dirname(__file__))}/results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])

adv_mx_sse = nc.groups['m_X_error']
adv_mx_sse = adv_mx_sse['total'][:]


# Conservative result:
con_dirname = f'{abspath(dirname(__file__))}/results/NL_sphere_conservative_order_{order}_consistency_ncells_{Ne}{extra_name}/field_output.nc'
conservative_nc = Dataset(f'{abspath(dirname(__file__))}/results/{con_dirname}/diagnostics.nc')

con_mx_sse = conservative_nc.groups['m_X_error']
con_mx_sse = con_mx_sse['total'][:]

# Plot these tracer densities
plt.figure()
plt.semilogy(times, adv_mx_sse, label='advective')
plt.semilogy(times, adv_mx_sse, label='conservative')
plt.xlabel('Time (s)')
plt.ylabel('|T_d(t) - T_d(0)| / T_d(0)')
plt.legend()

print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()