'''
DEV NOTE: This functions are deprecated and are thus not called by the primary integrator. These are 
mainly useful for testing on small meshes for internal consistency 


Simulation of fractional Brownian motion subject to a harmonic potential 
using standard reaction schemes

Several integrators are available, including Euler, RK4, and a fancy method that uses LSODA and 
which exploits the linearity of the diffusion operator for the system.


Notes:

This solver restricts itself to the case of absolute radial symmetry (the lowest-order spherical 
	harmonic) because of the symmetry of the harmonic potential being used

The solution to the radial equation includes a factor of 4 \pi r^2 and so it represents a radial
probability density of finding a random walker at a givne radius, rather than a concentration profile


William Gilpin, 2015

'''

from matplotlib.pyplot import *
from scipy import *
from numpy import *

from diffusion_integrator_funcs import *


def rk4_mesh(yinit, ynxt, settings, params):
    """
    ynxt : function
        The function that computes the next step in the forward march
    
    dt : double
        The spatial step
        
    dx : double
        The time step size
        
    use_crank : boolean
        Whether to use Crank-Nicolson forward marching, in which the next timestep
        is determined by equal weights of the previous two timesteps. If not, the
        integrator uses standard Forward Euler
        
    settings : dict
        A dictionary containing all the necessary steps and 
        objects for the integration
        
    params : dict
        A dictionary containing all the necessary parameters for the 
        differential equation system.
        
    """
    use_crank = True
    
    space = settings['space']
    times = settings['times']
    dt = times[2] - times[1]
    
    solmesh = zeros([len(yinit), len(times)])
    solmesh[:,0] = yinit
    tL = len(solmesh[0,:])
    for tind, time in enumerate(times):
        if tind > 0 and tind < tL:
            
            if use_crank == False:
                ind_list = [1]
                wgt = 1.0
            else:
                ind_list = [1,2]
                wgt = .5
            
            solmesh[:,tind] = solmesh[:,tind-1]
            for ii in ind_list:
                k1 = ynxt( solmesh[:,tind-ii], time, dt, settings, params)          
                k2 = ynxt( solmesh[:,tind-ii]+ 0.5*k1, time + 0.5*dt, dt, settings, params)
                k3 = ynxt( solmesh[:,tind-ii]+ 0.5*k2, time + 0.5*dt, dt, settings, params)
                k4 = ynxt( solmesh[:,tind-ii]+ k3, time+dt, dt, settings, params)
                solmesh[:,tind] += wgt*dt*(k1 + 2.0*(k2 + k3) + k4)/6.0
    
    # get gaussian part
    gauss = (1/sqrt(2*pi*params['POT_DIAM']**2))*exp(-space**2/(4*params['POT_DIAM']**2))
    gausspart = tile(gauss,(tL,1)).T
    solmesh = solmesh*gausspart
    
    # now convert to pdf
    radpart = tile(4*pi*(space**2),(tL,1)).T
    solmesh = radpart*solmesh
    return solmesh
    
def euler_mesh(yinit, ynxt, settings, params):
    """
    ynxt : function
        The function that returns the discrete derivative at a given timepoint
    
    dt : double
        The spatial step
        
    dx : double
        The time step size

    Notes: For almost any situation in which you might call this function, you might
    as well just use rk4_mesh. This code is slow and mainly exists for testing purposes
    """
    
    space = settings['space']
    times = settings['times']
    dt = times[2] - times[1]

    solmesh = zeros([len(yinit), len(times)])
    solmesh[:,0] = yinit
    tL = len(solmesh[0,:])
    for tind, time in enumerate(times):
        if tind > 0 and tind < tL:
            solmesh[:,tind] = solmesh[:,tind-1] + dt*ynxt( solmesh[:,tind-1], time, dt, settings, params)
    
    # get gaussian part
    # gauss = (1/sqrt(2*pi*params['POT_DIAM']**2))*exp(-space**2/(4*params['POT_DIAM']**2))
    # gausspart = tile(gauss,(tL,1)).T
    # solmesh = solmesh*gausspart
    
    # # now convert to pdf
    # radpart = tile(4*pi*(space**2),(tL,1)).T
    # solmesh = radpart*solmesh

    gauss = eq_dist(space, params['POT_DIAM'])
    gauss = gauss[:, None]
    sol = sol*gauss

    return solmesh





