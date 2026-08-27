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

from firedrake import Constant, ExtrudedMesh, as_vector, min_value, exp, cos
from gusto import *
import numpy as np

##########################################

# Specify the type of limiter to use
#limiter_type = 'none'
#limiter_type = 'standard'
limiter_type = 'mmr'

# IC type
#ic = 'Gaussian'
#ic = 'square'
ic = 'cylinder'

##########################################
# Parameters
ncells_per_edge = 24
radius = 6371220.      # radius of the sphere, in m
u_max = 10.       # Transporting velocity
dt = 450.    # Timestep size
tmax = 12*24*60*60.         # Simulation length
order = 1    # Next-to-lowest order elements
theta_c = 0.
lamda_c = 0.
dumpfreq = 20.
g_max = 0.5            # Maximum amplitude of Gaussian density perturbations
b0 = 2                 # Controls the width of the chemical blobs
rho_b = 1.

mesh = GeneralCubedSphereMesh(radius, ncells_per_edge, degree=2)
xyz = SpatialCoordinate(mesh)

# Only use order 1 elements
domain = Domain(mesh, dt, 'RTCF', 1)

# get lat lon coordinates
lamda, theta, _ = lonlatr_from_xyz(xyz[0], xyz[1], xyz[2])


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

#tracers = [m_X, rho_d]
tracers = [rho_d, m_X]

V = domain.spaces("HDiv")
eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu = V)

if limiter_type == 'none':
    dirname = f'mmr_test_no_limiter_sphere_{ic}'
elif limiter_type == 'standard':
    dirname = f'mmr_test_standard_limiter_sphere_{ic}'
elif limiter_type == 'mmr':
    dirname = f'mmr_test_mmr_limiter_sphere_{ic}'

# I/O
output = OutputParameters(
    dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=False
)

td = TracerDensity('m_X', 'rho_d', method='solve')

diagnostic_fields = [td, SteadyStateError('m_X')]

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
transport_methods = [DGUpwind(eqn, 'rho_d'), DGUpwind(eqn, 'm_X')]

#time_varying_velocity = False

time_varying_velocity=True
tau = tmax 
def u_t(t):
    k = 5.*radius/tau
    u_background = 2*pi*radius/tau
    lamda_prime = lamda - u_background*t

    u_zonal = (
        u_background*cos(theta)
        - k*(sin(lamda_prime/2)**2)*sin(2*theta)*(cos(theta)**2)*cos(pi*t/tau)
    )
    u_merid = 0.5*k*sin(lamda_prime)*(cos(theta)**3)*cos(pi*t/tau)

    return xyz_vector_from_lonlatr(u_zonal, u_merid, Constant(0.0), xyz)

stepper = PrescribedTransport(
    eqn, transport_scheme, io, time_varying_velocity, transport_methods
)

stepper.setup_prescribed_expr(u_t)

rho_d_0 = rho_b + 0.5*cos(theta)

X, Y, Z = xyz
X1, Y1, Z1 = xyz_from_lonlatr(lamda_c, theta_c, radius)

# Gaussian perturbation
if ic == 'Gaussian':
    m_X_0 = g_max*exp(-(b0/(radius**2))*((X-X1)**2 + (Y-Y1)**2 + (Z-Z1)**2))
elif ic == 'square':
# Square discontinuous shape - easier to get negatives
    m_X_0 = conditional((theta > -np.pi/10),
                        conditional(theta < np.pi/10,
                                    conditional(lamda > -np.pi/10,
                                                conditional(lamda < np.pi/10,
                                                            1.0,
                                                            0.0), 0.0), 0.0), 0.0)
elif ic == 'cylinder':
# Move to slotted cylinders!
    theta_c1 = 0.0         # latitude of first cylinder, in rad
    theta_c2 = 0.0         # latitude of second cylinder, in rad
    lamda_c1 = -pi/4       # longitude of first cylinder, in rad
    lamda_c2 = pi/4        # longitude of second cylinder, in rad

    m_X_0 = conditional(
                great_arc_angle(lamda, theta, lamda_c1, theta_c1) < 0.5,
                conditional(
                    abs(lamda - lamda_c1) < 1./12.,
                    conditional(theta - theta_c1 < -5./24., 1.0, 0.0),
                    1.0
                ),
                conditional(
                    great_arc_angle(lamda, theta, lamda_c2, theta_c2) < 0.5,
                    conditional(
                        abs(lamda - lamda_c2) < 1./12.,
                        conditional(theta - theta_c2 > 5./24., 1.0, 0.0),
                        1.0
                    ),
                    0.0
                )
            )


stepper.fields("m_X").interpolate(m_X_0)
stepper.fields("rho_d").interpolate(rho_d_0)

#u0 = stepper.fields("u")

# Solid body rotation wind of u_lon = u_0 cos(theta), u_lat = 0.
#u_max = Constant(2*pi*radius/(12*24*60*60))
#print(u_max)
#u0.project(as_vector([-u_max*Y/radius, u_max*X/radius, Constant(0.0)]))

# What about a non-divergent wind field from Nair, Lauritzen


stepper.run(t=0, tmax=tmax)