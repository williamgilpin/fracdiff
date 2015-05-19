'''

A set of helper functions for the fraction diffusion integrator code

William Gilpin, 2015

'''

from matplotlib.pyplot import *
from scipy import *
from numpy import *

def smoothstep(x,center=1.0,sharpness=1.0):
    '''
    return a smoothed step function that interpolates between one and zero
    '''
    return 1-1/(1+exp(-x+center))**sharpness



def lap2d(npts):
    """
    Build a central-difference laplace operator. Use a higher-order forward difference
    operator at the edges
    http://en.wikipedia.org/wiki/Finite_difference_coefficient
    """
    op = zeros([npts, npts])
    a = -2*ones((1, npts))[0]
    b = ones((1, npts-1))[0]
    op = diag(a, 0) + np.diag(b, -1) + np.diag(b, 1)
    op[0,1] = 2
    
    op[0,0] = 2
    op[0,1] = 5
    op[0,2] = 4
    op[0,3] = 1
    
    op[-1,-2] = 2
    
    return op


# @jit
def grad1D(npts):
    """
    Build a central-difference derivative operator. Use a higher-order forward difference
    operator at the edges
    http://en.wikipedia.org/wiki/Finite_difference_coefficient
    """
    op = zeros([npts, npts])
    a = zeros((1, npts))[0]
    b = .5*ones((1, npts-1))[0]
    op = diag(a, 0) + np.diag(b, -1) + np.diag(b, 1)
    
    op[0,0] = 0.
    op[0,1] = 1.
    
    op[0,0] = -1.5
    op[0,1] = 2.0
    op[0,2] = -.5
    
    op[-1,-1] = 1.
    op[-1,-2] = 0.
    return op


def make_fht(times, sol):
    '''
    Make an approximate first-passage time distribution given the survival probability distribution
    '''
    fht = (-diff(sum(sol, axis = 0)))
    return fht


def janky_laplace(signal, num_lvls = 50):
    '''
    finds the relative weights of various timescales present in
    a monotonically-decaying signal
    '''
    signal = signal.copy()
    signal = signal/ sum(signal)

    npts = len(signal)
    
#     # estimate average timescale present
#     krate = mean( diff(log(signal)) )
#     krate = -krate
    
# #     print ('mean timescale: ' + str(1/krate))
    
    xpts = arange(npts)
    
    min_timescale = 1
    max_timescale = npts
    
    kscales = 1/linspace(min_timescale, max_timescale, num_lvls)
     
    basis_vectors = exp(-outer(kscales, xpts))
    basis_vectors = (basis_vectors.T/sum(basis_vectors,axis=1)).T

    
    return kscales, basis_vectors


def expspace(lowlim, uplum, npts):
    '''
    logarithmically sample an interval
    '''
    return exp(linspace(log(lowlim), log(uplum),npts))

def gram_schmidt(X, row_vecs=True, norm = True):
    
    if not row_vecs:
        X = X.T
        
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
        
    if row_vecs:
        return Y
    else:
        return Y.T