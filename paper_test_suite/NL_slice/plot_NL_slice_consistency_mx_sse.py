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
#results_dir = '/data/home/ta440/firedrake_07_03_24/src/gusto/gusto_conservative_transport/paper_test_suite/NL_slice/'

# Or, for current directory:
results_dir = ''

# For any additional components to the name:
#extra_name = 'td_solve'
#extra_name = 'proj_18aug'
extra_name = ''


# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #

if order == 0:
    dxz = 200
    #dxz = 100
elif order == 1:
    dxz = 100
    #dxz = 50

plot_name = f'{abspath(dirname(__file__))}/figures/NL_slice_{extra_name}ord_{order}_consistency_m_l2_sse_over_time_dxz{dxz}.png'

# Advective result:
adv_dirname = f'NL_slice_{extra_name}advective_order_{order}_consistency_dxz_{dxz}'
nc = Dataset(f'{abspath(dirname(__file__))}/results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])

adv_mx_sse = nc.groups['m_X_error']
adv_mx_sse = adv_mx_sse['l2'][:]

# Normalise:
adv_mX = nc.groups['m_X']
adv_mX = adv_mX['l2'][:]
adv_mX0 = adv_mX[0]
adv_mx_sse = adv_mx_sse/adv_mX0


# Conservative result:
con_dirname = f'NL_slice_{extra_name}conservative_order_{order}_consistency_dxz_{dxz}'
conservative_nc = Dataset(f'{abspath(dirname(__file__))}/results/{con_dirname}/diagnostics.nc')

con_mx_sse = conservative_nc.groups['m_X_error']
con_mx_sse = con_mx_sse['l2'][:]

# Normalise:
con_mX = nc.groups['m_X']
con_mX = con_mX['l2'][:]
con_mX0 = con_mX[0]
con_mx_sse = con_mx_sse/con_mX0


font_opts = {'size': 14, 'family': 'serif'}
plt.rc('font', **font_opts)


# Plot these tracer densities
plt.figure()
plt.semilogy(times, adv_mx_sse, label='advective', c='r')
plt.semilogy(times, con_mx_sse, label='conservative', c='b')
plt.xlabel('Time (s)', size=16)
plt.ylabel('Relative mixing ratio error', size=16)
plt.xlim([0,2000])
plt.xticks([0,500, 1000, 1500, 2000])
plt.legend(loc='lower center', prop={'size': 16}, bbox_to_anchor=(0.5, -0.4))
plt.ylim([1e-16,1e-5])

plt.grid()

print(f'Saving figure to {plot_name}')
plt.savefig(plot_name, bbox_inches='tight')
plt.close()