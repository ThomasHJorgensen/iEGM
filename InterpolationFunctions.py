import numpy as np
from scipy import interpolate

def regressions_poly(nodes,degree):
    # returns nodes.size-by-degree+1 array of polynomials
    if not hasattr(nodes,'__len__'):
        nodes = np.array([nodes])
    num_nodes = nodes.size

    poly = np.ones((num_nodes,degree+1))
    for p in range(1,degree+1): # constant already first element
        poly[:,p] = nodes * poly[:,p-1]

    return poly

def regression_coefs(nodes,outcomes,degree):
    if not hasattr(nodes,'__len__'):
        nodes = np.array([nodes])
    num_nodes = nodes.size

    assert num_nodes>degree+1, f'Number of nodes ({num_nodes}) must be greater than polynomial degree ({degree}+1)!'

    # calculate polynomial
    poly = regressions_poly(nodes,degree)

    # regression coefficients
    return np.linalg.solve( poly.T @ poly , poly.T @ outcomes )


def regression_interp(nodes,coefs):
    if not hasattr(nodes,'__len__'):
        nodes = np.array([nodes])
    num_nodes = nodes.size

    # calculate polynomial
    degree = coefs.size-1
    poly = regressions_poly(nodes,degree)

    # predict outcome
    return poly @ coefs

def setup_Bspline(nodes,outcomes,par,degree=3):
    par.interp_Bspline = interpolate.splrep(nodes,outcomes,s=0,k=degree)

def interp_Bspline(x,par):
    return interpolate.BSpline(*par.interp_Bspline)(x)