from gusto import *
from firedrake import PeriodicIntervalMesh, ExtrudedMesh, Constant, ge, le, exp, cos, \
    sin, conditional, interpolate, SpatialCoordinate, VectorFunctionSpace, \
    Function, assemble, dx, FunctionSpace, pi, min_value, acos, as_vector

import numpy as np

u"""
Conservative Transport in Gusto:
Test Case 1 
  
This script implements a the test given in the 'Charney-Phillips trilemma'
paper by Bendall, Wood, Thuburn, and Cotter. This considers a planar version of
the Gaussian test case given by Nair and Laurtizen. 

This tests a coupled transport equation for moisture.

The mixing ratio obeys an advective transport equation:
∂/∂t (m_X) + (u.∇)m_X = 0 (1)

Whereas the dry density obeys the conservative form:
∂/∂t (ρ_d) + ∇.(ρ_d*u) = 0 (2)

In this script we test the ability to transport both of these variables
in a conservative manner. 

This means solving the following transport equation in place of (1),
∂/∂t (ρ_d*m_X) + ∇.(ρ_d*m_X*u) = 0 (3)

There are two configurations that can be run:
  The 'convergence' configuration has an initial condition of a linearly 
  varying density field and two Gaussian bumps for the mixing ratio.
  The 'consistency' configuration has an initial condition of a 
  constant mixing ratio and two Gaussian bumps for the density.

This test will have the mixing ratio and dry density in the
same function space, DG.

"""

# Specify whether to run the 'convergence' or 'consistency' version of the test.
case = 'convergence'

# Domain
Lx = 2000.
Hz = 2000.

# Time parameters
dt = 2.
tmax = 2000.

nlayers = 100.  # horizontal layers
columns = 100.  # number of columns

dx = Lx/nlayers
dz = Hz/columns

period_mesh = PeriodicIntervalMesh(columns, Lx)
mesh = ExtrudedMesh(period_mesh, layers=nlayers, layer_height=Hz/nlayers)
domain = Domain(mesh, dt, "CG", 1)
x,z = SpatialCoordinate(mesh)

# Choose spaces for the tracers
# Use DG to collocate the density and mixing ratio
rho_d_space = 'DG'
m_X_space = 'DG'

V_rho = domain.spaces(rho_d_space)
V_m_X = domain.spaces(m_X_space)

# Specify whether or not to use conservative form for the tracers.
conservative = True

# Define the mixing ratio and density as tracers

if conservative:
    m_X = ActiveTracer(name='m_X', space=m_X_space,
                     variable_type=TracerVariableType.mixing_ratio,
                     transport_eqn=TransportEquationType.tracer_conservative,
                     density_name='rho_d')
else:
    m_X = ActiveTracer(name='m_X', space=m_X_space,
                     variable_type=TracerVariableType.mixing_ratio,
                     transport_eqn=TransportEquationType.advective)
                 
rho_d = ActiveTracer(name='rho_d', space=rho_d_space,
                 variable_type=TracerVariableType.density,
                 transport_eqn=TransportEquationType.conservative)

tracers = [m_X,rho_d]

# Equation
V = domain.spaces("HDiv")

eqn = CoupledTransportEquation(domain, active_tracers=tracers, Vu = V)

if conservative:
    dirname = "test_1_fre"+case
else:
    dirname = "test_1_not_conservative_proj_"+case

# Dump the solution at each day
dumpfreq = int(100./dt)

# Set dump_nc = True to use tomplot.
output = OutputParameters(dirname=dirname,
                          dumpfreq = dumpfreq,
                          dump_nc = True,
                          dump_vtus = False)

# Use a tracer density diagnostic to track conservation
diagnostic_fields = [TracerDensity('m_X','rho_d',method='project')]

#diagnostic_fields = [assemble(rho_d*m_X*dx)]

io = IO(domain, output, diagnostic_fields=diagnostic_fields)

# Set up the divergent, time-varying, velocity field
U = Lx/tmax
W = U/10.

def u_t(t):
  xd = x - U*t
  u = U - (W*pi*Lx/Hz)*cos(pi*t/tmax)*cos(2*pi*xd/Lx)*cos(pi*z/Hz)
  w = 2*pi*W*cos(pi*t/tmax)*sin(2*pi*xd/Lx)*sin(pi*z/Hz)
  
  u_expr = as_vector((u,w))
  
  return u_expr

# Specify locations of the two Gaussians
xc1 = 5.*Lx/8.
zc1 = Hz/2.

xc2 = 3.*Lx/8.
zc2 = Hz/2.

def l2_dist(xc,zc):
  return min_value(abs(x-xc), Lx-abs(x-xc))**2 + (z-zc)**2

lc = 2.*Lx/25.
m0 = 0.02

# Set the initial state from the configuration choice
if case == 'convergence':
  f0 = 0.05
  
  rho_t = 0.5
  rho_b = 1.
  
  rho_d_0 = rho_b + z*(rho_t-rho_b)/Hz  
  
  g1 = f0*exp(-l2_dist(xc1,zc1)/(lc**2))
  g2 = f0*exp(-l2_dist(xc2,zc2)/(lc**2))
  
  m_X_0 = m0 + g1 + g2
  
elif case == 'consistency':
  f0 = 0.5
  rho_b = 0.5
  
  g1 = f0*exp(-l2_dist(xc1,zc1)/(lc**2))
  g2 = f0*exp(-l2_dist(xc2,zc2)/(lc**2))
  
  rho_d_0 = rho_b + g1 + g2
  
  m_X_0 = m0 + 0*x
  
else:
  raise NotImplementedError('Specified case is not recognised.')

# This is just for the 'reference' solution
# with advective form transport
apply_limiter = False

if conservative:
    # Use the mass-weighted evaluations for timestepping.
    transport_scheme = SSPRK3(domain, increment_form=False)
elif apply_limiter:
    sublimiters = {'m_X': DG1Limiter(V_m_X), 
                   'rho_d': DG1Limiter(V_rho)}
    MixedLimiter = MixedFSLimiter(eqn, sublimiters)

    transport_scheme = SSPRK3(domain, limiter=MixedLimiter)
else:
    transport_scheme = SSPRK3(domain)

transport_methods = [DGUpwind(eqn, "m_X"), DGUpwind(eqn, "rho_d")]
    
time_varying = True
    
# Timestepper
stepper = SplitPrescribedTransport(eqn, transport_scheme, 
                                   io, time_varying, transport_methods
)

#stepper = SplitPrescribedTransport(eqn, transport_scheme, io, prescribed_transporting_velocity=True, spatial_methods=transport_methods)

stepper.setup_prescribed_expr(u_t)

# Initial Conditions
stepper.fields("m_X").interpolate(m_X_0)
stepper.fields("rho_d").interpolate(rho_d_0)
u0 = stepper.fields("u")
u0.project(u_t(0))

stepper.run(t=0, tmax=tmax)
