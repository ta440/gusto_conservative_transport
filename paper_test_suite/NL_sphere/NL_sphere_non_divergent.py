u"""
A test for the Conservative Transport in Gusto paper (Tim and Tom):

'NL_sphere'.

This implements a test from the Nair Laurtizen paper:
'A class of deformational flow test cases for linear transport
problems on the sphere'. 

We will test consistency and conservation with co-located function
spaces, where rho and m both lie in DG.
A cubed-sphere mesh is used so that we have quadrilateral elements.

- The 'convergence' configuration has an initial condition of a linearly
      varying density field and two Gaussian bumps for the mixing ratio.
- The 'consistency' configuration has an initial condition of a
      constant mixing ratio and two Gaussian bumps for the density.

We will test this setup with order 0 and order 1 finite elements.

This script prescribes a divergence-free velocity field,
which is defined through the streamfunction.

"""


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from firedrake import (
    exp, cos, sin, SpatialCoordinate, Constant, pi
)
from gusto import *

NL_sphere_defaults = {
    'conservative_transport': False,  # whether to use conservative transport
    'configuration': 'convergence',   # 'convergence or 'consistency'
    'order': 1,                       # order of the finite element spaces
    'ncells_per_edge': 24,            # num points per cubed sphere panel edge
    'dt': 450.0,                      # 7.5 minutes, which should be sufficient for ref level 5.
    'tmax': 12.*24.*60.*60.,          # 12 days
    'dumpfreq': 288,                  # 8 outputs: once every 1.5 days
}


def NL_sphere(
        conservative_transport=NL_sphere_defaults['conservative_transport'],
        configuration=NL_sphere_defaults['configuration'],
        order=NL_sphere_defaults['order'],
        ncells_per_edge=NL_sphere_defaults['ncells_per_edge'],
        dt=NL_sphere_defaults['dt'],
        tmax=NL_sphere_defaults['tmax'],
        dumpfreq=NL_sphere_defaults['dumpfreq']
):

    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #

    tau = 12.*24.*60.*60.  # time period for reversible flow, in s
    radius = 6371220.      # radius of sphere, in m
    theta_c1 = 0.0         # latitude of first blob, in rad
    theta_c2 = 0.0         # latitude of second blob, in rad
    lamda_c1 = -pi/4       # longitude of first blob, in rad
    lamda_c2 = pi/4        # longitude of second blob, in rad
    rho_b = 1              # Base dry density
    m0 = 0.02              # Base mixing ratio value

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    print('Using conservative transport?: ', conservative_transport)

    # Domain
    mesh = GeneralCubedSphereMesh(radius, ncells_per_edge, degree=2)
    xyz = SpatialCoordinate(mesh)
    domain = Domain(mesh, dt, 'RTCF', order)

    # Use DG for both tracers (colocated)
    tracer_space = 'DG'

    # Define the mixing ratio and density as tracers
    rho_d = ActiveTracer(
        name='rho_d', space=tracer_space,
        variable_type=TracerVariableType.density,
        transport_eqn=TransportEquationType.conservative
    )

    if conservative_transport:
        m_X = ActiveTracer(
            name='m_X', space=tracer_space,
            variable_type=TracerVariableType.mixing_ratio,
            transport_eqn=TransportEquationType.tracer_conservative,
            density_name='rho_d'
        )
    else:
        m_X = ActiveTracer(
            name='m_X', space=tracer_space,
            variable_type=TracerVariableType.mixing_ratio,
            transport_eqn=TransportEquationType.advective
        )

    tracers = [rho_d, m_X]

    # Equation
    V = domain.spaces("HDiv")
    eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu=V)

    if conservative_transport:
        transport_type = 'conservative'
    else:
        transport_type = 'advective'

    dirname = f'NL_sphere_{transport_type}_order_{str(order)}_{configuration}_ncells_{str(ncells_per_edge)}_dt_{dt}_non_div_flow_proj_low'

    # I/O
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=False
    )

    # Use a tracer density diagnostic to track conservation
    # Compare interpolate vs solve methods.
    td = TracerDensity('m_X', 'rho_d')
    #td = TracerDensity('m_X', 'rho_d', method='solve')

    diagnostic_fields = [
        td, SteadyStateError('m_X'),
        ZonalComponent('u'), MeridionalComponent('u')
    ]

    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    # Details of transport
    transport_methods = [DGUpwind(eqn, "m_X"), DGUpwind(eqn, "rho_d")]

    if order == 1:
        suboptions = {}
    elif order == 0:
        VCG1 = FunctionSpace(mesh, 'CG', 1)
        #VDG1 = domain.spaces('DG1_equispaced')    
        VDG1 = domain.spaces('DG1_equispaced')    
    
        if conservative_transport:
            suboptions = {'rho_d': RecoveryOptions(embedding_space=VDG1,
                                                   recovered_space=VCG1,
                                                   project_low_method='recover'),
                          'm_X': ConservativeRecoveryOptions(embedding_space=VDG1,
                                                             recovered_space=VCG1,
                                                             project_low_method='conservative_project',
                                                             project_high_method='conservative_project',
                                                             rho_name='rho_d',
                                                             orig_rho_space=domain.spaces(tracer_space))
                                                             }
        else:
            suboptions = {'rho_d': RecoveryOptions(embedding_space=VDG1,
                                                   recovered_space=VCG1,
                                                   project_low_method='project'),
                          'm_X': RecoveryOptions(embedding_space=VDG1,
                                                 recovered_space=VCG1,
                                                 project_low_method='project')
                        }
    else:
        raise NotImplementedError('Higher-order spaces are not'
                                  + 'implemented for this test case.')

    opts = MixedFSOptions(suboptions=suboptions)

    if conservative_transport:
        transport_scheme = SSPRK3(domain, options=opts, rk_formulation=RungeKuttaFormulation.predictor)
    else:
        transport_scheme = SSPRK3(domain, options=opts)

    # Time stepper
    time_varying_velocity = True
    stepper = PrescribedTransport(
        eqn, transport_scheme, io, time_varying_velocity, transport_methods
    )

    # ------------------------------------------------------------------------ #
    # Initial conditions
    # ------------------------------------------------------------------------ #

    # Transporting wind ------------------------------------------------------ #
    lamda, theta, _ = lonlatr_from_xyz(xyz[0], xyz[1], xyz[2])

    H1 = domain.spaces('H1')
    psi = Function(H1)
    u0 = stepper.fields("u")

    k = 10.*radius/tau
    lamda_prime = lamda - 2*pi*stepper.t/tau

    # Divergence-free wind, obtained from stream function
    psi_expr = radius*(
        k*((sin(lamda_prime)*cos(theta))**2)*cos(pi*stepper.t/tau)
        - 2.*pi*radius*sin(theta)/tau
    )

    u_expr = domain.perp(grad(psi))

    psi_interpolator = Interpolator(psi_expr, psi)
    u_projector = Projector(u_expr, u0)

    # Set up the non-divergent, time-varying, velocity field
    def apply_prescribed_velocity(t):
        psi_interpolator.interpolate()
        u_projector.project()
        return

    stepper.setup_prescribed_apply(apply_prescribed_velocity)

    # Density and mixing ratio ------------------------------------------------ #

    # Specify locations of the two bumps
    X = cos(theta)*cos(lamda)
    Y = cos(theta)*sin(lamda)
    Z = sin(theta)

    X1 = cos(theta_c1)*cos(lamda_c1)
    Y1 = cos(theta_c1)*sin(lamda_c1)
    Z1 = sin(theta_c1)

    X2 = cos(theta_c2)*cos(lamda_c2)
    Y2 = cos(theta_c2)*sin(lamda_c2)
    Z2 = sin(theta_c2)

    # Define the two Gaussian bumps
    g1 = exp(-5*((X-X1)**2 + (Y-Y1)**2 + (Z-Z1)**2))
    g2 = exp(-5*((X-X2)**2 + (Y-Y2)**2 + (Z-Z2)**2))

    if configuration == 'convergence':
        g_max = 0.05
        rho_d_0 = rho_b + 0.5*cos(theta)
        m_X_0 = m0 + g_max*g1 + g_max*g2
    elif configuration == 'consistency':
        g_max = 0.5
        rho_d_0 = rho_b + g_max*g1 + g_max*g2
        m_X_0 = m0 + 0*cos(theta)
    else:
        raise ValueError('Specified configuration is not valid')

    # Set fields
    stepper.fields("m_X").interpolate(m_X_0)
    stepper.fields("rho_d").interpolate(rho_d_0)

    # ------------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------------ #

    stepper.run(t=0, tmax=tmax)


# ---------------------------------------------------------------------------- #
# MAIN
# ---------------------------------------------------------------------------- #


if __name__ == "__main__":

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--conservative_transport',
        help="Whether to apply conservative transport or not",
        type=bool,
        default=NL_sphere_defaults['conservative_transport']
    )
    parser.add_argument(
        '--configuration',
        help="The test configuration to use, 'convergence' or 'consistency'",
        type=str,
        default=NL_sphere_defaults['configuration']
    )
    parser.add_argument(
        '--order',
        help="The order of the finite elements",
        type=int,
        default=NL_sphere_defaults['order']
    )
    parser.add_argument(
        '--ncells_per_edge',
        help="The number of cells per edge of the icosahedron",
        type=int,
        default=NL_sphere_defaults['ncells_per_edge']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=NL_sphere_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=NL_sphere_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=NL_sphere_defaults['dumpfreq']
    )
    args, unknown = parser.parse_known_args()

    NL_sphere(**vars(args))
