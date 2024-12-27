
u"""
Conservative Transport in Gusto:

This script runs the Laurtizen et al. (2015) Terminator Toy
test case. This examines the interaction of two species
in the transport equation. There is coupling 
between the two species to model combination 
and dissociation.

There is a dry density, rho_d, and two mixing ratios,
X and X2. We want to test using conservative form
for both of these mixing ratios. Specifically, the mixing ratios
will live in the theta space, with rho_d in DG.

We will only use order 1 elements in this test,
but will examine conservation when the tracers are in 
colocated and staggered spaces.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from firedrake import IcosahedralSphereMesh, Constant, ge, le, exp, cos, \
    sin, conditional, interpolate, SpatialCoordinate, VectorFunctionSpace, \
    Function, assemble, dx, FunctionSpace, pi, max_value, acos, as_vector

from gusto import *
import numpy as np

terminator_toy_defaults = {
    'conservative_transport': False, # Whether to use conservative transport
    'diff_spaces': False,            # Colocated or staggered spaces
    'ncells_per_edge': 24,           # num points per cubed sphere panel edge
    'dt': 450.0,                     # 7.5 minutes
    'tmax': 12.*24.*60.*60.,         # 12 days
    'dumpfreq': 288,                 # once every 1.5 days with default values
}

def terminator_toy(
        conservative_transport=terminator_toy_defaults['conservative_transport'],
        diff_spaces=terminator_toy_defaults['diff_spaces'],
        ncells_per_edge=terminator_toy_defaults['ncells_per_edge'],
        dt=terminator_toy_defaults['dt'],
        tmax=terminator_toy_defaults['tmax'],
        dumpfreq=terminator_toy_defaults['dumpfreq']
):


    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #

    tau = tmax             # time period of reversible wind, in s
    radius = 6371220.      # radius of the sphere, in m
    theta_cr = pi/9.       # central latitude of first reaction rate, in rad
    lamda_cr = -pi/3.      # central longitude of first reaction rate, in rad
    k1_max = 1.            # amplitude of first reaction rate parameter, in 1/s
    k2 = 1.                # second reaction rate parameter, in 1/s
    theta_c1 = 0.          # central latitude of first chemical blob, in rad
    theta_c2 = 0.          # central latitude of second chemical blob, in rad
    lamda_c1 = -pi/4.      # central longitude of first chemical blob, in rad
    lamda_c2 = pi/4.       # central longitude of second chemical blob, in rad
    b0 = 5                 # controls the width of the chemical blobs

    # Domain
    print('Using conservative transport?: ', conservative_transport)

    print('Using different spaces for rho and X,X2? ', diff_spaces)

    # Domain
    mesh = GeneralCubedSphereMesh(radius, ncells_per_edge, degree=2)
    xyz = SpatialCoordinate(mesh)
    domain = Domain(mesh, dt, 'RTCF', 1)

    # get lat lon coordinates
    lamda, theta, _ = lonlatr_from_xyz(xyz[0], xyz[1], xyz[2])

    # Use staggered spaces
    rho_d_space = 'DG'

    if diff_spaces:
        m_X_space = 'theta'
    else:
        m_X_space = 'DG'

    V_rho = domain.spaces(rho_d_space)
    V_m_X = domain.spaces(m_X_space)


    # Define the dry density and the two species as tracers
    rho_d = ActiveTracer(name='rho_d', space=rho_d_space,
                    variable_type=TracerVariableType.density,
                    transport_eqn=TransportEquationType.conservative)

    if conservative_transport:
        X_tracer = ActiveTracer(name='X_tracer', space=m_X_space,
                        variable_type=TracerVariableType.mixing_ratio,
                        transport_eqn=TransportEquationType.tracer_conservative,
                        density_name='rho_d')
        
        X2_tracer = ActiveTracer(name='X2_tracer', space=m_X_space,
                        variable_type=TracerVariableType.mixing_ratio,
                        transport_eqn=TransportEquationType.tracer_conservative,
                        density_name='rho_d')
    else:
        X_tracer = ActiveTracer(name='X_tracer', space=m_X_space,
                        variable_type=TracerVariableType.mixing_ratio,
                        transport_eqn=TransportEquationType.advective)
        
        X2_tracer = ActiveTracer(name='X2_tracer', space=m_X_space,
                        variable_type=TracerVariableType.mixing_ratio,
                        transport_eqn=TransportEquationType.advective)

    tracers = [rho_d, X_tracer, X2_tracer]

    # Equation
    V = domain.spaces("HDiv")

    eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu = V)

    if conservative_transport:
        transport_type='conservative'
    else:
        transport_type='advective'

    if diff_spaces:
        dirname = 'terminator_toy_diff_spaces_'+transport_type+'_ncells_'+str(ncells_per_edge)
    else:
        dirname = 'terminator_toy_same_spaces_'+transport_type+'_ncells_'+str(ncells_per_edge)
   


    # Set dump_nc = True to use tomplot.
    output = OutputParameters(dirname=dirname,
                              dumpfreq = dumpfreq,
                              dump_nc = True,
                              dump_vtus = False)     

    X_mass = TracerDensity('X_tracer', 'rho_d')
    X2_mass = TracerDensity('X2_tracer', 'rho_d')

    io = IO(domain, output, diagnostic_fields = [X_mass, X2_mass])

    k1 = max_value(0, sin(theta)*sin(theta_cr) + cos(theta)*cos(theta_cr)*cos(lamda-lamda_cr))
    k2 = 1

    terminator_stepper = BackwardEuler(domain)

    physics_schemes = [(TerminatorToy(eqn, k1=k1, k2=k2, species1_name='X_tracer',
                        species2_name='X2_tracer'), terminator_stepper)]

    X, Y, Z = xyz
    X1, Y1, Z1 = xyz_from_lonlatr(lamda_c1, theta_c1, radius)
    X2, Y2, Z2 = xyz_from_lonlatr(lamda_c2, theta_c2, radius)

    g1 = exp(-(b0/(radius**2))*((X-X1)**2 + (Y-Y1)**2 + (Z-Z1)**2))
    g2 = exp(-(b0/(radius**2))*((X-X2)**2 + (Y-Y2)**2 + (Z-Z2)**2))

    rho_expr = g1 + g2

    X_T_0 = 4e-6
    r = k1/(4*k2)
    D_val = sqrt(r**2 + 2*X_T_0*r)

    # Initial condition for each species
    X_0 = D_val - r
    X2_0 = 0.5*(X_T_0 - D_val + r)

    def u_t(t):
        k = 10*radius/tau

        u_zonal = (
            k*(sin(lamda - 2*pi*t/tau)**2)*sin(2*theta)*cos(pi*t/tau)
            + ((2*pi*radius)/tau)*cos(theta)
        )
        u_merid = k*sin(2*(lamda - 2*pi*t/tau))*cos(theta)*cos(pi*t/tau)

        return xyz_vector_from_lonlatr(u_zonal, u_merid, Constant(0.0), xyz)

    if diff_spaces:
        Vt_brok = FunctionSpace(mesh, BrokenElement(V_m_X.ufl_element()))
        
        if conservative_transport:
            suboptions = {'rho_d': EmbeddedDGOptions(embedding_space=Vt_brok),
                          'X_tracer':ConservativeEmbeddedDGOptions(project_back_method='conservative_project',
                                                            rho_name='rho_d',
                                                            orig_rho_space=V_rho),
                          'X2_tracer':ConservativeEmbeddedDGOptions(project_back_method='conservative_project',
                                                             rho_name='rho_d',
                                                             orig_rho_space=V_rho)
                                                             }
        else:
            suboptions = {'rho_d':RecoveryOptions(embedding_space=Vt_brok,
                                                  recovered_space=V_m_X,
                                                  project_low_method='recover'),
                          'X_tracer': EmbeddedDGOptions(),
                          'X2_tracer': EmbeddedDGOptions()}
    else:
        suboptions = {}

    opts = MixedFSOptions(suboptions=suboptions)

    if conservative_transport:
        transport_scheme = SSPRK3(domain, options=opts, rk_formulation=RungeKuttaFormulation.predictor)
    else:
        transport_scheme = SSPRK3(domain, options=opts)

    transport_method = [DGUpwind(eqn, 'rho_d'), DGUpwind(eqn, 'X_tracer'), DGUpwind(eqn, 'X2_tracer')]
                                        
    time_varying_velocity=True
    stepper = SplitPrescribedTransport(eqn, transport_scheme, io, 
                                        time_varying_velocity,
                                        spatial_methods=transport_method,
                                        physics_schemes = physics_schemes)
                                        
    stepper.setup_prescribed_expr(u_t)

    # Initial conditions
    stepper.fields("rho_d").interpolate(rho_expr)
    stepper.fields("X_tracer").interpolate(X_0)
    stepper.fields("X2_tracer").interpolate(X2_0)

    # Run until Termination!
    stepper.run(t=0, tmax=tmax)


if __name__ == "__main__":

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--conservative_transport',
        help="Whether to apply conservative transport or not",
        type=bool,
        default=terminator_toy_defaults['conservative_transport']
    )
    parser.add_argument(
        '--diff_spaces',
        help="Whether to use different spaces for the tracers or not",
        type=bool,
        default=terminator_toy_defaults['diff_spaces']
    )
    parser.add_argument(
        '--ncells_per_edge',
        help="The number of cells per edge of the icosahedron",
        type=int,
        default=terminator_toy_defaults['ncells_per_edge']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=terminator_toy_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=terminator_toy_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=terminator_toy_defaults['dumpfreq']
    )
    args, unknown = parser.parse_known_args()

    terminator_toy(**vars(args))
