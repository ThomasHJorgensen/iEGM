import numpy as np

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




def Chebyshev_nodes(points,num_nodes):
    # end-points
    a = points[0]
    b = points[-1]

    # calculate the Chebyshev nodes
    nodes = np.nan + np.zeros((num_nodes))  # Initial vector to store the Chebyshev nodes
    for k in range(num_nodes):                
        
        # Step1: Compute the m Chebyshev interpolation notes in [-1,1]    
        zk = Chebyshev_interp_notes(k,num_nodes)

        # Step 2: Adjust the nodes to the [a,b] interval
        nodes[k] = (zk+1)*((b-a)/2)+a

    return nodes


def Chebyshev_coefs(f_actual,num_nodes,degree):

    assert (num_nodes>=degree+1), 'The specified parameters are not acceptable. Make sure m>n'

    coefs = np.nan + np.ones(degree+1)
    for i in range(degree+1):             # Loop over the degree of the approximation polynomial. 
        nom = 0                           # Initial value for step 4
        denom = 0                         # Initial value for step 4
        for k in range(num_nodes):                # Loop over the approximation notes
            
            # Step1: Compute the m Chebyshev interpolation notes in [-1,1]    
            zk = Chebyshev_interp_notes(k,num_nodes)

            # Step 2: Compute Chebyshev coefficients. Tn=cos(i*acos(zk))
            nom += f_actual[k]*Tn(zk,i)
            denom += Tn(zk,i)**2
            if k==num_nodes-1:
                coefs[i] = nom/denom

    return coefs

def interp_chebyshev(points,coefs,degree,x_eval):

    a = points[0]
    b = points[-1]

    f_approx = 0.0                             # Initial value of the approximating function
    for i in range(degree+1):                  # Loop over the degree of the approximation polynomial. 
        f_approx += coefs[i]*Tn(2*(x_eval-a)/(b-a)-1,i)  # The Chebyshev approximation of f(x)

    return f_approx


def Chebyshev_illustration(fhandle,points,m,n):
    
    # This is the Chebyshev Interpolation (Regression algorithm)      
    #  in approximation of a scalar function, f(x):R->R                
    #    The approach follow Judd (1998, Allgortihm 6.2, p. 223)         
#############################################################################
# INPUT ARGUMENTS:
#             fhandle:               The funtion, that should be approximated
#             interval:              The interval for the approximation of f(x).
#             m:                     number of nodes used to construct the approximation. NOTE: m>=n+1
#             n:                     Degree of approximation-polynomial
# 
# OUTPUT ARGUMENTS:
#             f_approx:              The vector of approximated function values
#             f_actual:              The vector of actual function values
#             points:                The vector of points, for which the function is approximated
##################################################################################

    assert (m>=n+1), 'The specified parameters are not acceptable. Make sure m>n'

    a = points[0]
    b = points[-1]
    number = points.size
    f_approx = np.nan + np.zeros((number))  # Initial vector to store the approximated function values
    f_actual = np.nan + np.zeros((number))  # Initial vector to store the actual function values

    for x in range(number):                   # Loop over the x values
        ai = np.nan +np.zeros((n+1))         # Initial vector to store the Chebyshev coefficients
        f_hat = 0                             # Initial value of the approximating function
        for i in range(n+1):                  # Loop over the degree of the approximation polynomial. 
            nom = 0                           # Initial value for step 4
            denom = 0                         # Initial value for step 4
            for k in range(m):                # Loop over the approximation notes
                
                # Step1: Compute the m Chebyshev interpolation notes in [-1,1]    
                zk = -np.cos(((2*(k+1)-1)/(2*m))*np.pi)

                # Step 2: Adjust the nodes to the [a,b] interval
                xk = (zk+1)*((b-a)/2)+a

                # Step 3: Evaluate f at the approximation notes.
                yk = fhandle(xk);  

                # Step 4: Compute Chebyshev coefficients. Tn=cos(i*acos(zk))
                nom += yk*Tn(zk,i)
                denom += Tn(zk,i)**2
                if k==m-1:
                    ai[i] = nom/denom
            
            f_hat = f_hat+ai[i]*Tn(2*(points[x]-a)/(b-a)-1,i)  # The Chebyshev approximation of f(x)
            f_temp = fhandle(points[x])                       # Actual function value, f(x)

        f_approx[x] = f_hat
        f_actual[x] = f_temp

    return f_approx, f_actual
def Chebyshev_interp_notes(k,num_nodes):
    return -np.cos(((2*(k+1)-1)/(2*num_nodes))*np.pi)

def Tn(x,n):
    return np.cos(n*np.arccos(x))


