"""
Conservative Transport in Gusto:
Test Case 4

The moist rising bubble test from Bryan & Fritsch (2002), in a cloudy
atmosphere.

The rise of the thermal is fueled by latent heating from condensation.
These use the compressible Euler equations. We additionally define
transport equations for water vapour (w_v) and cloud water (c_w).
When these are defined conservatively, we have the following five equations:

∂u/∂t + (u.∇)u + 2Ω×u + c_p*θ*∇Π + g = 0,                                 
∂ρ/∂t + ∇.(ρ*u) = 0,                                                      
∂θ/∂t + (u.∇)θ = 0,                                                       
∂(ρ*w_v)/∂t + ∇.(ρ*w_v*u) = 0, 
∂(ρ*c_w)/∂t + ∇.(ρ*c_w*u) = 0.

where Π is the Exner pressure, g is the gravitational vector, Ω is the
planet's rotation vector and c_p is the heat capacity of dry air at constant

"""

from gusto import *
from gusto import thermodynamics
from firedrake import (PeriodicIntervalMesh, ExtrudedMesh,
                       SpatialCoordinate, conditional, cos, pi, sqrt,
                       NonlinearVariationalProblem,
                       NonlinearVariationalSolver, TestFunction, dx,
                       TrialFunction, Function,
                       LinearVariationalProblem, LinearVariationalSolver)
import sys

# ---------------------------------------------------------------------------- #
# Test case parameters
# ---------------------------------------------------------------------------- #

dt = 1.0
L = 10000.
H = 10000.

deltax = 200
tmax = 1000.

# ---------------------------------------------------------------------------- #
# Set up model objects
# ---------------------------------------------------------------------------- #

# Domain
nlayers = int(H/deltax)
ncolumns = int(L/deltax)

m = PeriodicIntervalMesh(ncolumns, L)
mesh = ExtrudedMesh(m, layers=nlayers, layer_height=H/nlayers)
degree = 1
domain = Domain(mesh, dt, 'CG', degree)

# Equation
params = CompressibleParameters()

# Choose whether to use conservative transport
conservative=True

# By default, water and clour vapour are in theta.
# Rho is in L2.
if conservative:
    tracers = [WaterVapour(transport_eqn=TransportEquationType.tracer_conservative,
                           density_name='rho'),
               CloudWater(transport_eqn=TransportEquationType.tracer_conservative,
                          density_name='rho')]
else:
    tracers = [WaterVapour(), CloudWater()]

eqns = CompressibleEulerEquations(domain, params, active_tracers=tracers)

# I/O
if conservative:
    dirname = 'test_4_conservative_mixed_opts'
else:
    dirname = 'test_4_not_conservative'

dumpfreq = int(tmax/(10.*dt))

# Set dump_nc = True to use tomplot.
output = OutputParameters(dirname=dirname,
                          dumpfreq = dumpfreq,
                          dump_nc = True,
                          dump_vtus = False)

diagnostic_fields = [Theta_e(eqns), \
                     TracerDensity('water_vapour', 'rho'),
                     TracerDensity('cloud_water', 'rho')]

io = IO(domain, output, diagnostic_fields=diagnostic_fields)

# Transport schemes
# Use Recovery for rho to have it embedded in the same spaces, also?
# rho is in L2, u is in H(div)

# Use increment_form for any explicit multistage timestepping
# for conservative transport of tracers.

                                                   
suboptions = {'water_vapour': EmbeddedDGOptions(),
              'cloud_water': EmbeddedDGOptions()}
mixed_opts = MixedFSOptions(suboptions=suboptions)
                          
if conservative:
    #transported_fields = [SSPRK3(domain, ["rho", "water_vapour", "cloud_water"], increment_form=False),     
    transported_fields = [SSPRK3(domain, ["rho", "water_vapour", "cloud_water"], options=mixed_opts, rk_formulation=RungeKuttaFormulation.predictor),
                          SSPRK3(domain, "theta", options=EmbeddedDGOptions()),
                          TrapeziumRule(domain, "u")]                         
                          
                          
else:
    transported_fields = [SSPRK3(domain, "rho"),
                          SSPRK3(domain, "theta", options=EmbeddedDGOptions()),
                          SSPRK3(domain, "water_vapour", options=EmbeddedDGOptions()),
                          SSPRK3(domain, "cloud_water", options=EmbeddedDGOptions()),
                          TrapeziumRule(domain, "u")]

transport_methods = [DGUpwind(eqns, field) for field in ["u", "rho", "theta", "water_vapour", "cloud_water"]]

# Linear solver
linear_solver = CompressibleSolver(eqns)

# Physics schemes (condensation/evaporation)
physics_schemes = [(SaturationAdjustment(eqns), ForwardEuler(domain))]

# Time stepper
stepper = SemiImplicitQuasiNewton(eqns, io, transported_fields,
                                  transport_methods,
                                  linear_solver=linear_solver,
                                  physics_schemes=physics_schemes)

# ---------------------------------------------------------------------------- #
# Initial conditions
# ---------------------------------------------------------------------------- #

u0 = stepper.fields("u")
rho0 = stepper.fields("rho")
theta0 = stepper.fields("theta")
water_v0 = stepper.fields("water_vapour")
water_c0 = stepper.fields("cloud_water")

# spaces
Vu = domain.spaces("HDiv")
Vt = domain.spaces("theta")
Vr = domain.spaces("DG")
x, z = SpatialCoordinate(mesh)
quadrature_degree = (4, 4)
dxp = dx(degree=(quadrature_degree))

# Define constant theta_e and water_t
Tsurf = 320.0
total_water = 0.02
theta_e = Function(Vt).assign(Tsurf)
water_t = Function(Vt).assign(total_water)

# Calculate hydrostatic fields
saturated_hydrostatic_balance(eqns, stepper.fields, theta_e, water_t)

# make mean fields
theta_b = Function(Vt).assign(theta0)
rho_b = Function(Vr).assign(rho0)
water_vb = Function(Vt).assign(water_v0)
water_cb = Function(Vt).assign(water_t - water_vb)
exner_b = thermodynamics.exner_pressure(eqns.parameters, rho_b, theta_b)
Tb = thermodynamics.T(eqns.parameters, theta_b, exner_b, r_v=water_vb)

# define perturbation
xc = L / 2
zc = 2000.
rc = 2000.
Tdash = 2.0
r = sqrt((x - xc) ** 2 + (z - zc) ** 2)
theta_pert = Function(Vt).interpolate(
    conditional(r > rc,
                0.0,
                Tdash * (cos(pi * r / (2.0 * rc))) ** 2))

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

stepper.set_reference_profiles([('rho', rho_b),
                                ('theta', theta_b),
                                ('water_vapour', water_vb),
                                ('cloud_water', water_cb)])

# ---------------------------------------------------------------------------- #
# Run
# ---------------------------------------------------------------------------- #

stepper.run(t=0, tmax=tmax)