'''
Test the mean mixing ratio limiter.
Do so using a transport equation.
Although 1d would be cheaper, 2d is 
easier to implement at the moment.


Test three different configurations:
1) No limiting
2) Standard DG1 limiter that is not mass conserving
3) MMR limiter

Idea for setup is a constant density,
and have a step function for the mixing ratio. This can
easily generate negatives from errors in the transport
scheme.

Consider a density and mixing ratio both in the DG space.

'''

from firedrake import PeriodicIntervalMesh, Constant, ExtrudedMesh, as_vector, min_value, exp
from gusto import *
import numpy as np

##########################################

# Specify the type of limiter to use
#limiter_type = 'none'
#limiter_type = 'standard'
limiter_type = 'mmr'

# test_type = 'Gaussian'
#test_type = 'square'

##########################################
# Parameters
ncells_1d = 20
Lx = 100      # Length of the domain
initial_wind = 10.      # Transporting velocity
dt = 0.05    # Timestep size
tmax = 10.       # Simulation length
order = 1    # Next-to-lowest order elements
xc = Lx/2
zc = Lx/2
lc = 2.*Lx/25.   # Decay rate of Gaussian
dumpfreq = 20.
f0 = 0.0        # Base tracer value
f_pert = 0.1

period_mesh = PeriodicIntervalMesh(ncells_1d, Lx)
mesh = ExtrudedMesh(period_mesh, layers=ncells_1d, layer_height=Lx/ncells_1d)
domain = Domain(mesh, dt, "CG", order)
x, z = SpatialCoordinate(mesh)

tracer_space = 'DG'
V_tracer = domain.spaces(tracer_space)

rho_d = ActiveTracer(name='rho_d', space=tracer_space,
                   variable_type=TracerVariableType.density,
                   transport_eqn=TransportEquationType.conservative)

m_X = ActiveTracer(
            name='m_X', space=tracer_space,
            variable_type=TracerVariableType.mixing_ratio,
            transport_eqn=TransportEquationType.tracer_conservative,
            density_name='rho_d'
        )

tracers = [rho_d, m_X]

V = domain.spaces("HDiv")
eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu = V)

if limiter_type == 'none':
    dirname = f'mmr_test_no_limiter_vary_rho'
elif limiter_type == 'standard':
    dirname = f'mmr_test_standard_limiter_vary_rho'
elif limiter_type == 'mmr':
    dirname = f'mmr_test_mmr_limiter_vary_rho_subtract_mean'

# I/O
output = OutputParameters(
    dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=False
)

td_X = TracerDensity('m_X', 'rho_d', method='solve')

diagnostic_fields = [td_X, SteadyStateError('m_X')]

io = IO(domain, output, diagnostic_fields=diagnostic_fields)


solver_parameters = conservative_tracer_parameters(V_tracer, num_fields=2)

if limiter_type == 'none':
    transport_scheme = SSPRK3(domain, rk_formulation=RungeKuttaFormulation.predictor,
                              solver_parameters=solver_parameters)
elif limiter_type == 'standard':
    sublimiters = {'m_X': DG1Limiter(V_tracer), 
                   'rho_d': DG1Limiter(V_tracer)}
    MixedLimiter = MixedFSLimiter(eqn, sublimiters)
    transport_scheme = SSPRK3(domain, rk_formulation=RungeKuttaFormulation.predictor,
                              limiter=MixedLimiter, solver_parameters=solver_parameters)
elif limiter_type == 'mmr':
    augmentation = MeanMixingRatio(domain, eqn, ['m_X'])
    transport_scheme = SSPRK3(domain, augmentation=augmentation,
                              rk_formulation=RungeKuttaFormulation.predictor,
                              solver_parameters=solver_parameters)

# Details of transport
transport_methods = [DGUpwind(eqn, "rho_d"), DGUpwind(eqn, "m_X")]

time_varying_velocity = False
stepper = PrescribedTransport(
    eqn, transport_scheme, io, time_varying_velocity, transport_methods
)
#stepper = SplitPrescribedTransport(
#    eqn, transport_scheme, io, time_varying_velocity, transport_methods, physics_schemes=None
#)

#u_t = u0
#stepper.setup_prescribed_expr(u_t)

#rho_d_0 = Constant(0.5)

# Linearly varying density field
rho_d_0 = Constant(0.5) + cos(z*pi/Lx)**2

#m_X_0 = conditional((x > 40.),
#                     conditional(x < 60.,
#                                 conditional(z > 40.,
#                                             conditional(z < 60.,
#                                                         1.0,
#                                                         0), 0), 0), 0)

def l2_dist(xc,zc):
    return min_value(abs(x-xc), Lx-abs(x-xc))**2 + (z-zc)**2

m_X_0 = f0 + f_pert*exp(-l2_dist(xc, zc)/lc**2)

stepper.fields("m_X").interpolate(m_X_0)
stepper.fields("rho_d").interpolate(rho_d_0)
u0 = stepper.fields("u")
u0.project(as_vector([initial_wind, 0]))

stepper.run(t=0, tmax=tmax)