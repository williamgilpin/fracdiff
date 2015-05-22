'''
Core operators and functions for the fractional Brownian Integrator using LSODA

This file includes a series of nested functions that is eventually passed as a single matrix,
the jacobian, to the low-level version of the scipy odeint function, which allows the integration
method and jacobian treatment to be custom-selected depending on the type of problem.

William Gilpin, 2015
'''

from matplotlib.pyplot import *
from scipy import *
from numpy import *

from diffusion_integrator_funcs import *

def eq_dist(pts, pot_width):
    '''
    The radial gaussian distribution for my parametricization of the diffusion equation.
    Compare the (4a^2) term to the equilibrium distribution in Doi and Edwards in order
    to find a in terms of constants like N and the Kuhn length 
    '''
    
    gausspart = exp(-(pts**2)/(4*pot_width**2))
    radpart = 4*pi*pts**2
    norm_fac = 8 * ((pot_width)**3) * (pi**(3/2))
    
    return gausspart*radpart/norm_fac


def nxt_step(timepoint, yvals, settings_and_params):
    """
    Return a space-discretized approximation of the differential equation
    at the given timepoint. Our equation happens to be linear, and so the next
    step is given by a simple dot product with the relevant operator
    """
    settings, params = settings_and_params
    jac = radiff_timeop(timepoint, yvals, [settings, params])
    nxt_vals = jac.dot(yvals)
    return nxt_vals


def time_nondim(time, params):
    '''
    Given a span of time values or a timestep, return the rescaled value
    of time used for our nondimensional form of the reaction-diffusion equation
    
    params : dict
        An array of parameter values for the STANDARD, UNIT-BEARING diffusion equation
        
    time : array
        A double array of time values (or just a single double)
    
    '''
    
    DCOEFF = params['DCOEFF']
    ALPHA = params['ALPHA']

    
    return DCOEFF*(time**ALPHA)



def radiff_spaceop(yvals, space, pot_diam):
    '''
    The diffusion operator in a harmonic potential.
    
    Drop the time evolution operator from a radial diffusion equation in a 
    harmonic potential and just calculate results of the spatial operator 
    on a given snapshot of the concentration profile
    
    Note: The diffusion coefficient is not included in the operator, that part
    goes in the time evolution function
    
    space : array
        The x values on which the simulation occurs
    
    yvals : the current concentration profile on the space
    
    pot_diam : the one parameter of the diffusion operator
    
    William Gilpin, 2015
    '''
    
    a = pot_diam
    allspace = space
    dx = allspace[2]-allspace[1]
    
    L = len(space)
    
    # second derivative step
    lap = lap2d(L)/dx**2
    jac = lap
    
    
    # first derivative steps
    drv = grad1D(L)/dx
    jac += (4*drv/allspace) - (3/(4*a**2))*allspace*drv
    
    
    # other steps
    iden = identity(L)
    jac += (1/(8*a**4))*(allspace**2)*iden
    jac += -(1/a**2)*iden
    # assume \ell = 0
    jac += (2/(allspace**2))*iden
    
    return jac



def radiff_timeop(timepoint, yvals, settings_and_params):
    '''
    Time evolution operator for three-dimensional, fractional brownian
    motion in a harmonci potential with a reactive well.
    
    '''
    
    settings, params = settings_and_params
    
    allspace = settings['space']
    
    WELL_DIAM = params['WELL_DIAM']
    POT_DIAM = params['POT_DIAM']
    KAPPA = params['KAPPA']
    DCOEFF = params['DCOEFF']
    ALPHA = params['ALPHA']
    
    DCOEFF = ALPHA*(timepoint**(ALPHA-1))*DCOEFF
    
    # remove this and change allspace
    a = POT_DIAM

    L = len(allspace)
    iden = identity(L)
    
    jac = radiff_spaceop(yvals, allspace, a)
    # now multiply by diffusion coefficient and all that
    jac = -DCOEFF*jac
    
    # now implement reaction step
    drain_window = double(allspace < WELL_DIAM)
    
    allinds = linspace(0,len(allspace),len(allspace))
    drain_window = smoothstep(allinds,center = sum(drain_window),sharpness=3)

    jac += -KAPPA*drain_window*iden
    
    return jac



def radiff_timeop_nondim(timepoint, yvals, settings_and_params):
    """
    Time operator form of a reaction-fractional diffusion equation
    where the parameters have been nondimensionalized. This function
    is similar in design to its dimensional cousin radiff_timeop
    """
    
    settings, params = settings_and_params
    
    allspace = settings['space']
    dx = allspace[2]-allspace[1]
    
    WELL_DIAM = params['WELL_DIAM']
    POT_DIAM = params['POT_DIAM']
    KAPPA = params['KAPPA']
    ALPHA = params['ALPHA']
    DCOEFF = params['DCOEFF']
    
    # Rescale into dimensionless form
    KAPPA = ( KAPPA/( ALPHA*(DCOEFF**(1/ALPHA)) ) )
    
    
    # remove this and change allspace
    a = POT_DIAM

    L = len(allspace)
    iden = identity(L)
    
    jac = radiff_spaceop(yvals, allspace, a)
    
    # now multiply by diffusion coefficient and all that
    jac = -jac
    
    # now implement reaction step
    drain_window = double(allspace < WELL_DIAM)
    
    allinds = linspace(0,len(allspace),len(allspace))
    drain_window = smoothstep(allinds,center = sum(drain_window),sharpness=3)

    jac += -KAPPA*drain_window*iden
    
    return jac



# def jacob(timepoint, yvals, settings_and_params):
#     """
#     Operator form of a system of coupled differential equations in one variable
#     """
    
#     settings, params = settings_and_params
    
#     allspace = settings['space']
#     dx = allspace[2]-allspace[1]
    
#     WELL_DIAM = params['WELL_DIAM']
#     POT_DIAM = params['POT_DIAM']
#     KAPPA = params['KAPPA']
#     DCOEFF = params['DCOEFF']
#     ALPHA = params['ALPHA']
    
#     DCOEFF = ALPHA*(timepoint**(ALPHA-1))*DCOEFF
    
#     # remove this and change allspace
#     a = POT_DIAM

#     L = len(space)
    
#     # second derivative step
#     lap = lap2d(L)/dx**2
#     jac = lap
    
    
#     # first derivative steps
#     drv = grad1D(L)/dx
#     jac += (4*drv/allspace) - (3/(4*a**2))*allspace*drv
    
    
#     # other steps
#     iden = identity(L)
#     jac += (1/(8*a**4))*(allspace**2)*iden
#     jac += -(1/a**2)*iden
#     # assume \ell = 0
#     jac += (2/(allspace**2))*iden
    
    
#     # now multiply by diffusion coefficient and all that
#     jac = -DCOEFF*jac
    
#     # now implement reaction step
#     drain_window = double(allspace < WELL_DIAM)
    
#     allinds = linspace(0,len(allspace),len(allspace))
#     drain_window = smoothstep(allinds,center = sum(drain_window),sharpness=3)

#     jac += -KAPPA*drain_window*iden
    
    
#     return jac


# def nxt_step(yvals, timepoint, dt, settings, params):
#     """
#     Return a space-discretized approximation of the differential equation
#     at the given timepoint. Our equation happens to be linear, and so the next
#     step is given by a simple dot product with the relevant operator
#     """
#     jac = radiff_timeop(timepoint, yvals, [settings, params])
#     nxt_vals = jac.dot(yvals)
#     return nxt_vals
