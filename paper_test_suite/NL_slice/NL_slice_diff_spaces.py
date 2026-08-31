u"""
A test for the Conservative Transport in Gusto paper (Tim and Tom):

'NL_slice'.
  
This implements the test in the 'Charney-Phillips trilemma'
paper by Bendall, Wood, Thuburn, and Cotter, which is a planar version of
the Gaussian test case given by Nair and Laurtizen. 

We will test consistency and conservation with staggered function 
spaces, where rho is in DG and m is in the theta space.

- The 'convergence' configuration has an initial condition of a linearly
      varying density field and two Gaussian bumps for the mixing ratio.
- The 'consistency' configuration has an initial condition of a
      constant mixing ratio and two Gaussian bumps for the density.

We will test this setup with order 0 and order 1 finite elements.

"""


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from firedrake import (
    PeriodicIntervalMesh, ExtrudedMesh, cos, exp, sin, SpatialCoordinate, pi,
    min_value, as_vector, FunctionSpace, BrokenElement
)
from gusto import *


NL_slice_defaults = {
    'conservative_transport': False,   # whether to use conservative transport
    'configuration': 'convergence',   # 'convergence or 'consistency'
    'order': 1,                       # order of the finite element spaces                
    'ncells_1d': 100,                 # number of points in x and z directions
    'dt': 2.0,                        # Fixed dt for all tests
    'tmax': 2000.,
    'dumpfreq': 125,                  # with defaults gives eight outputs
}

def NL_slice(
        conservative_transport=NL_slice_defaults['conservative_transport'],
        configuration=NL_slice_defaults['configuration'],
        order=NL_slice_defaults['order'],
        ncells_1d=NL_slice_defaults['ncells_1d'],
        dt=NL_slice_defaults['dt'],
        tmax=NL_slice_defaults['tmax'],
        dumpfreq=NL_slice_defaults['dumpfreq']
):


    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #

    tau = tmax       # time period of reversible flow, in s
    Lx = 2000.       # width of domain, in m
    Hz = 2000.       # height of domain, in m
    w_factor = 0.1   # deformational factor, dimensionless
    xc1 = 5.*Lx/8.   # x-coordinate of centre of Gaussian bump 1
    zc1 = Hz/2.      # z-coordinate of centre of Gaussian bump 1
    xc2 = 3.*Lx/8.   # x-coordinate of centre of Gaussian bump 2
    zc2 = Hz/2.      # z-coordinate of centre of Gaussian bump 2
    lc = 2.*Lx/25.   # Decay rate of Gaussian
    m0 = 0.02        # Base mixing ratio value

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    print('Using conservative transport?: ', conservative_transport)

    period_mesh = PeriodicIntervalMesh(ncells_1d, Lx)
    mesh = ExtrudedMesh(period_mesh, layers=ncells_1d, layer_height=Hz/ncells_1d)
    domain = Domain(mesh, dt, "CG", order)
    x, z = SpatialCoordinate(mesh)

    # Use staggered spaces for the tracers
    rho_d_space = 'DG'
    m_X_space = 'theta'
    
    V_rho = domain.spaces(rho_d_space)
    V_m_X = domain.spaces(m_X_space)
    
    # Define the mixing ratio and density as tracers
    rho_d = ActiveTracer(
        name='rho_d', space=rho_d_space,
        variable_type=TracerVariableType.density,
        transport_eqn=TransportEquationType.conservative
    )
    
    if conservative_transport:
        m_X = ActiveTracer(
            name='m_X', space=m_X_space,
            variable_type=TracerVariableType.mixing_ratio,
            transport_eqn=TransportEquationType.tracer_conservative,
            density_name='rho_d'
        )
    else:
        m_X = ActiveTracer(
            name='m_X', space=m_X_space,
            variable_type=TracerVariableType.mixing_ratio,
            transport_eqn=TransportEquationType.advective
        )

    tracers = [rho_d, m_X]

    # Equation
    V = domain.spaces("HDiv")
    eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu = V)
    
    if conservative_transport:
        transport_type='conservative'
    else:
        transport_type='advective'
    
    dirname = 'NL_slice_diff_spaces_'+transport_type+'_order_'+str(order)+'_'+configuration+'_dxz_'+str(ncells_1d)
    
    # I/O
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=False
    )
    
    # Use a tracer density diagnostic to track conservation.
    # Use the solve method for the best accuracy
    td = TracerDensity('m_X', 'rho_d', method='solve')

    diagnostic_fields = [
        td, SteadyStateError('m_X'),
        XComponent('u'), ZComponent('u')
    ]

    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    # Details of transport
    transport_methods = [DGUpwind(eqn, "m_X"), DGUpwind(eqn, "rho_d")]
    
    # Specify options depending on the order of the space:
    if order == 0:
        # Specify recovery options for both tracers
        VCG1 = FunctionSpace(mesh, 'CG', 1)
        VDG1 = domain.spaces('DG1_equispaced')

        if conservative_transport:
            suboptions = {'rho_d': RecoveryOptions(embedding_space=VDG1,
                                                   recovered_space=VCG1,
                                                   project_low_method='project',
                                                   boundary_method=BoundaryMethod.taylor),
                          'm_X': ConservativeRecoveryOptions(embedding_space=VDG1,
                                                             recovered_space=VCG1,
                                                             boundary_method=BoundaryMethod.taylor,
                                                             project_low_method='conservative_project',
                                                             project_high_method='conservative_project',
                                                             rho_name='rho_d',
                                                             orig_rho_space=V_rho)
                                                             }
        else:
            suboptions = {'rho_d': RecoveryOptions(embedding_space=VDG1,
                                                   recovered_space=VCG1,
                                                   project_low_method='project',
                                                   boundary_method=BoundaryMethod.taylor),
                          'm_X': RecoveryOptions(embedding_space=VDG1,
                                                 recovered_space=VCG1,
                                                 project_low_method='project',
                                                 boundary_method=BoundaryMethod.taylor)
            }

    elif order != 1:
        raise NotImplementedError('Higher-order spaces have not been'
                                  + 'implemented for this test case.')

    if conservative_transport:

        transport_scheme = SSPRK3(domain, options=opts, rk_formulation=RungeKuttaFormulation.predictor)
    else:
        transport_scheme = SSPRK3(domain)
    
    
    time_varying_velocity = True
    stepper = PrescribedTransport(
        eqn, transport_scheme, io, time_varying_velocity, transport_methods
    )
    
    # Transporting wind ------------------------------------------------------ #
    # Set up the divergent, time-varying, velocity field
    def u_t(t):
        umax = Lx/tau
        xd = x - umax*t

        u = umax - (w_factor*umax*pi*Lx/Hz)*cos(pi*t/tau)*cos(2*pi*xd/Lx)*cos(pi*z/Hz)
        w = 2*pi*w_factor*umax*cos(pi*t/tau)*sin(2*pi*xd/Lx)*sin(pi*z/Hz)

        return as_vector([u, w])
    
    stepper.setup_prescribed_expr(u_t)
    
    # Initial condition
    
    def l2_dist(xc,zc):
      return min_value(abs(x-xc), Lx-abs(x-xc))**2 + (z-zc)**2
      
    if configuration == 'convergence':
        f0 = 0.05

        rho_t = 0.5
        rho_b = 1.

        rho_d_0 = rho_b + z*(rho_t-rho_b)/Hz

        g1 = f0*exp(-l2_dist(xc1, zc1)/lc**2)
        g2 = f0*exp(-l2_dist(xc2, zc2)/lc**2)

        m_X_0 = m0 + g1 + g2

    elif configuration == 'consistency':
        f0 = 0.5
        rho_b = 0.5

        g1 = f0*exp(-l2_dist(xc1, zc1)/lc**2)
        g2 = f0*exp(-l2_dist(xc2, zc2)/lc**2)

        rho_d_0 = rho_b + g1 + g2

        m_X_0 = m0 + 0*x

    else:
        raise ValueError('Specified configuration is not valid')
        
    # Set fields
    stepper.fields("m_X").interpolate(m_X_0)
    stepper.fields("rho_d").interpolate(rho_d_0)
    u0 = stepper.fields("u")
    u0.project(u_t(0))

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
        default=NL_slice_defaults['conservative_transport']
    )
    parser.add_argument(
        '--configuration',
        help="The test configuration to use, 'convergence' or 'consistency'",
        type=str,
        default=NL_slice_defaults['configuration']
    )
    parser.add_argument(
        '--order',
        help="The order of the finite elements",
        type=int,
        default=NL_slice_defaults['order']
    )
    parser.add_argument(
        '--ncells_1d',
        help="The number of cells in the x and y directions",
        type=int,
        default=NL_slice_defaults['ncells_1d']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=NL_slice_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=NL_slice_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=NL_slice_defaults['dumpfreq']
    )
    args, unknown = parser.parse_known_args()

    NL_slice(**vars(args))
