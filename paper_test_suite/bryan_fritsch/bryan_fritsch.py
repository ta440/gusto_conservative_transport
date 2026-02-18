"""

A test for the Conservative Transport in Gusto paper (Tim and Tom):

'bryan_fritsch'.

The moist rising bubble test from Bryan & Fritsch, 2002:
``A Benchmark Simulation for Moist Nonhydrostatic Numerical Models'', GMD.

The test simulates a rising thermal in a cloudy atmosphere, which is fueled by
latent heating from condensation.

This test will be run with order 0 and order 1 finite element spaces.

The two tracers are defined in the theta space 
and are transported conservatively.

In this version, we use a linearly varying mixing ratio distribution
to make it more difficult to ensure conservation.

"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from firedrake import (
    PeriodicIntervalMesh, ExtrudedMesh, SpatialCoordinate, conditional, cos, pi,
    sqrt, NonlinearVariationalProblem, NonlinearVariationalSolver, TestFunction,
    dx, TrialFunction, Function, as_vector, LinearVariationalProblem,
    LinearVariationalSolver, Constant, BrokenElement, assemble
)

from gusto import *

moist_bryan_fritsch_defaults = {
    'conservative_transport': False, # Whether to use conservative transport
    'order': 1, # Order of the finite elements
    'ncolumns': 50,
    'nlayers': 50,
    'dt': 2.0,
    'tmax': 1002.0,
    'dumpfreq': 50.
}


def moist_bryan_fritsch(
        conservative_transport=moist_bryan_fritsch_defaults['conservative_transport'],
        order=moist_bryan_fritsch_defaults['order'],
        ncolumns=moist_bryan_fritsch_defaults['ncolumns'],
        dt=moist_bryan_fritsch_defaults['dt'],
        tmax=moist_bryan_fritsch_defaults['tmax'],
        dumpfreq=moist_bryan_fritsch_defaults['dumpfreq']
):

    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #
    nlayers = ncolumns        # Set same resolution in x and z
    domain_width = 10000.     # domain width, in m
    domain_height = 10000.    # domain height, in m
    zc = 2000.                # vertical centre of bubble, in m
    rc = 2000.                # radius of bubble, in m
    Tdash = 2.0               # strength of temperature perturbation, in K
    Tsurf = 320.0             # background theta_e value, in K
    m0 = 0.02                 # Base mixing ratio, in kg/kg
    delta_m = 0.005           # Linear variation in mixing ratio with height

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    # Discretisation of advective term
    if order == 0:
        u_eqn_type = "vector_advection_form"
    else:
        u_eqn_type = "vector_invariant_form"

    # Domain
    base_mesh = PeriodicIntervalMesh(ncolumns, domain_width)
    mesh = ExtrudedMesh(
        base_mesh, layers=nlayers, layer_height=domain_height/nlayers
    )
    domain = Domain(mesh, dt, 'CG', order)


    # Set up the tracers and their transport schemes
    V_rho = domain.spaces('DG')
    V_theta = domain.spaces('theta')

    if conservative_transport:
        tracers = [WaterVapour(space='theta',
                               transport_eqn=TransportEquationType.tracer_conservative,
                               density_name='rho'),
                   CloudWater(space='theta',
                              transport_eqn=TransportEquationType.tracer_conservative,
                              density_name='rho')]
    else:
        tracers = [WaterVapour(), CloudWater()]

    # Equation
    params = CompressibleParameters(mesh)
    eqns = CompressibleEulerEquations(
        domain, params, active_tracers=tracers, u_transport_option=u_eqn_type
    )

    # I/O
    if conservative_transport:
        transport_type='conservative'
        print('Using conservative transport for the tracers')
    else:
        transport_type='advective'
        print('Not using conservative transport for the tracers')
    
    dirname = 'bryan_fritsch_feb18_2026_'+transport_type+'_order_'+str(order)+'dxz'+str(nlayers)
    
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_vtus=False, dump_nc=True
    )
    diagnostic_fields = [Theta_e(eqns), XComponent('u'), ZComponent('u'),
                         TracerDensity('water_vapour', 'rho'),
                         TracerDensity('cloud_water', 'rho'),
                         Sum('water_vapour','cloud_water'),
                         TracerDensity('water_vapour_plus_cloud_water', 'rho')]
    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    # Set up transport schemes
    if order == 0:
        VDG1 = domain.spaces("DG1_equispaced")
        VCG1 = FunctionSpace(mesh, "CG", 1)
        Vu_DG1 = VectorFunctionSpace(mesh, VDG1.ufl_element())
        Vu_CG1 = VectorFunctionSpace(mesh, "CG", 1)

        u_opts = RecoveryOptions(embedding_space=Vu_DG1,
                                 recovered_space=Vu_CG1,
                                 boundary_method=BoundaryMethod.taylor)
        theta_opts = RecoveryOptions(embedding_space=VDG1,
                                     recovered_space=VCG1,
                                     boundary_method=BoundaryMethod.taylor)

        if conservative_transport:
            suboptions = {'rho': RecoveryOptions(embedding_space=VDG1,
                                                recovered_space=VCG1,
                                                boundary_method=BoundaryMethod.taylor),
                         'water_vapour': ConservativeRecoveryOptions(embedding_space=VDG1,
                                                                     recovered_space=VCG1,
                                                                     rho_name="rho",
                                                                     orig_rho_space=V_rho,
                                                                     boundary_method=BoundaryMethod.taylor),
                         'cloud_water': ConservativeRecoveryOptions(embedding_space=VDG1,
                                                                    recovered_space=VCG1,
                                                                    rho_name="rho",
                                                                    orig_rho_space=V_rho,
                                                                    boundary_method=BoundaryMethod.taylor)}
        else:
            rho_opts = RecoveryOptions(embedding_space=VDG1,
                                       recovered_space=VCG1,
                                       boundary_method=BoundaryMethod.taylor)
            wv_opts = RecoveryOptions(embedding_space=VDG1,
                                      recovered_space=VCG1,
                                      boundary_method=BoundaryMethod.taylor)
            wc_opts = RecoveryOptions(embedding_space=VDG1,
                                      recovered_space=VCG1,
                                      boundary_method=BoundaryMethod.taylor)
    else:
        theta_opts = EmbeddedDGOptions()
        if conservative_transport:
            Vt_brok = FunctionSpace(mesh, BrokenElement(V_theta.ufl_element()))
            suboptions = {'rho': EmbeddedDGOptions(embedding_space=Vt_brok),
                          'water_vapour': ConservativeEmbeddedDGOptions(embedding_space=Vt_brok,
                                                                      rho_name="rho",
                                                                      orig_rho_space=V_rho),
                          'cloud_water': ConservativeEmbeddedDGOptions(embedding_space=Vt_brok,
                                                                     rho_name="rho",
                                                                     orig_rho_space=V_rho)}
        else:
            rho_opts = None
            wv_opts = EmbeddedDGOptions()
            wc_opts = EmbeddedDGOptions()

    transported_fields = [SSPRK3(domain, "theta", options=theta_opts)]
    
    if conservative_transport:
        mixed_opts = MixedFSOptions(suboptions=suboptions)
        transported_fields.append(SSPRK3(domain, ["rho", "water_vapour", "cloud_water"], options=mixed_opts, rk_formulation=RungeKuttaFormulation.predictor))
    else:
        transported_fields.append(SSPRK3(domain, 'rho', options=rho_opts))
        transported_fields.append(SSPRK3(domain, 'water_vapour', options=wv_opts))
        transported_fields.append(SSPRK3(domain, 'cloud_water', options=wc_opts))

    if order == 0:
        transported_fields.append(SSPRK3(domain, 'u', options=u_opts))
    else:
        transported_fields.append(TrapeziumRule(domain, 'u'))

    transport_methods = [
        DGUpwind(eqns, field) for field in
        ["u", "rho", "theta", "water_vapour", "cloud_water"]
    ]

    # Physics schemes (condensation/evaporation)
    physics_schemes = [(SaturationAdjustment(eqns), ForwardEuler(domain))]

    # Time stepper
    # Use 16 outer loops, 1 inner loop, to avoid
    # the linear solver breaking mass conservation
    stepper = SemiImplicitQuasiNewton(
        eqns, io, transported_fields, transport_methods,
        final_physics_schemes=physics_schemes, num_outer=16, num_inner=1
    )

    # ------------------------------------------------------------------------ #
    # Initial conditions
    # ------------------------------------------------------------------------ #

    u0 = stepper.fields("u")
    rho0 = stepper.fields("rho")
    theta0 = stepper.fields("theta")
    water_v0 = stepper.fields("water_vapour")
    water_c0 = stepper.fields("cloud_water")

    # spaces
    Vt = domain.spaces("theta")
    Vr = domain.spaces("DG")
    x, z = SpatialCoordinate(mesh)
    quadrature_degree = (4, 4)
    dxp = dx(degree=(quadrature_degree))

    # Define constant theta_e
    theta_e = Function(Vt).assign(Tsurf)

    # Initialise a linearly varying mixing_ratio
    total_water = m0 - z*delta_m/domain_height
    water_t = Function(Vt).interpolate(total_water)

    # Calculate hydrostatic fields
    saturated_hydrostatic_balance(eqns, stepper.fields, theta_e, water_t)

    # make mean fields
    theta_b = Function(Vt).assign(theta0)
    rho_b = Function(Vr).assign(rho0)
    water_vb = Function(Vt).assign(water_v0)
    water_cb = Function(Vt).assign(water_t - water_vb)

    # define perturbation
    xc = domain_width / 2
    r = sqrt((x - xc) ** 2 + (z - zc) ** 2)
    theta_pert = Function(Vt).interpolate(
        conditional(
            r > rc,
            0.0,
            Tdash * (cos(pi * r / (2.0 * rc))) ** 2
        )
    )

    # define initial theta
    theta0.interpolate(theta_b * (theta_pert / 300.0 + 1.0))

    # find perturbed rho
    gamma = TestFunction(Vr)
    rho_trial = TrialFunction(Vr)
    a = gamma * rho_trial * dxp
    L = gamma * (rho_b * theta_b / theta0) * dxp
    rho_problem = LinearVariationalProblem(a, L, rho0)
    rho_solver = LinearVariationalSolver(rho_problem)
    rho_solver.solve()

    # find perturbed water_v
    w_v = Function(Vt)
    phi = TestFunction(Vt)
    rho_averaged = Function(Vt)
    rho_recoverer = Recoverer(rho0, rho_averaged)
    rho_recoverer.project()

    exner = thermodynamics.exner_pressure(eqns.parameters, rho_averaged, theta0)
    p = thermodynamics.p(eqns.parameters, exner)
    T = thermodynamics.T(eqns.parameters, theta0, exner, r_v=w_v)
    w_sat = thermodynamics.r_sat(eqns.parameters, T, p)

    w_functional = (phi * w_v * dxp - phi * w_sat * dxp)
    w_problem = NonlinearVariationalProblem(w_functional, w_v)
    w_solver = NonlinearVariationalSolver(w_problem)
    w_solver.solve()

    water_v0.assign(w_v)
    water_c0.assign(water_t - water_v0)

    # wind initially zero
    u0.project(as_vector([Constant(0.0), Constant(0.0)]))

    stepper.set_reference_profiles(
        [
            ('rho', rho_b),
            ('theta', theta_b),
            ('water_vapour', water_vb),
            ('cloud_water', water_cb)
        ]
    )

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
        default=moist_bryan_fritsch_defaults['conservative_transport']
    )
    parser.add_argument(
        '--order',
        help="The order of the finite elements.",
        type=int,
        default=moist_bryan_fritsch_defaults['order']
    )
    parser.add_argument(
        '--ncolumns',
        help="The number of columns in the vertical slice mesh.",
        type=int,
        default=moist_bryan_fritsch_defaults['ncolumns']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=moist_bryan_fritsch_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=moist_bryan_fritsch_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=moist_bryan_fritsch_defaults['dumpfreq']
    )
    args, unknown = parser.parse_known_args()

    moist_bryan_fritsch(**vars(args))
