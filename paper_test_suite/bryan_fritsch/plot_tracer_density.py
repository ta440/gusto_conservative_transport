# A script to plot the tracer density from the 
# bryan fritsch bubble test

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
from tomplot import (set_tomplot_style, plot_convergence,
                     only_minmax_ticklabels, tomplot_legend_ax,
                     tomplot_legend_fig)
import numpy as np

order = 0

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/figures'
plot_name = f'{plot_dir}/bryan_frisch_ord_{order}_Td_over_time.png'


# Advective result:
adv_dirname = 'bryan_fritsch_advective_order_'+str(order)
nc = Dataset(f'results/{adv_dirname}/diagnostics.nc')

times = np.asarray(nc['time'])

T_wv = nc.groups['TracerDensity_water_vapour_rho']
Twv = T_wv['total'][:]
T_cw = nc.groups['TracerDensity_cloud_water_rho']
Tcw = T_cw['total'][:]

Td_adv = Twv + Tcw
Td_adv_diff = (Td_adv - Td_adv[0])/Td_adv[0]

# Conservative result:
con_dirname = 'bryan_fritsch_conservative_order_'+str(order)
conservative_nc = Dataset(f'results/{con_dirname}/diagnostics.nc')

T_wv = nc.groups['TracerDensity_water_vapour_rho']
Twv = T_wv['total'][:]
T_cw = nc.groups['TracerDensity_cloud_water_rho']
Tcw = T_cw['total'][:]

Td_con = Twv + Tcw
Td_con_diff = (Td_con - Td_con[0])/Td_con[0]

# Plot these tracer densities
plt.figure()
plt.plot(times, Td_adv_diff, label='advective')
plt.plot(times, Td_con_diff, label='conservative')
plt.xlabel('Time (s)')
plt.ylabel('|T_d(t) - T_d(0)| / T_d(0)')
plt.legend()