### Importing the libraries
import math
import Library.assign3 as a3


#Bracket function to find the interval
def Bracket(a,b,func,t,d):
    """Brackets the root of a function.

    Args:
        a (float): lower bound of the interval.
        b (float): upper bound of the interval.
        func (class): function for which we want to find the root.
        t ( int): maximum number of iterations.
        d (float): shifting parameter.
    Returns:
       Interval : [a,b] where the root lies.
    """
    x = func(a)
    y = func(b)
    if t == 10: 
        return 0
    if x*y < 0:
        print("a=",a,",b=",b,"\nIterations:",t,"\n")
        return a,b
    t+=1
    if x*y > 0:
        if abs(x) < abs(y):   
          return Bracket(float(a-d*(b-a)),b,func,t,d) 
        elif abs(x) > abs(y):
          return Bracket(a,float(b+d*(b-a)),func,t,d)
      
      
#### bisection method to find the root of a function

def bisection(f, a, b, tol=1e-6, maxiter=1000):
    """Find the root of a function using the bisection method.
    
    Parameters
    ----------
    f : function
        The function to find the root of.
    a : float
        The lower bound of the interval.
    b : float
        The upper bound of the interval.
    tol : float, optional
        The tolerance for the root.
    maxiter : int, optional
        The maximum number of iterations.
    
    Returns
    -------
    float
        The root of the function.
    """
    # create a list to store the values of x
    x = []    
    # check that the bounds are valid
    if f(a) * f(b) > 0:
        raise ValueError("The function must have opposite signs at the bounds.")
    
    # check that the bounds are in the correct order
    if a > b:
        a, b = b, a
    
    # iterate until the root is found
    for i in range(maxiter):
        # find the midpoint
        c = (a + b) / 2
        x.append(c)
        # check if the root is found
        if abs(f(c)) < tol:
            return x
        
        # check if the root is in the left interval
        if f(a) * f(c) < 0:
            b = c
        # check if the root is in the right interval
        else:
            a = c
    # if the root is not found, raise an error
    raise RuntimeError("The root was not found within the maximum number of iterations.")

#### Regula falsi method of root finding

def regula_falsi(f, a, b, tol=1e-6, max_iter=100):
    """Regula falsi method for finding roots of a function.
    Parameters
    ----------
    f : function
        The function for which we want to find a root.
    a : float
        Left endpoint of the interval.
    b : float
        Right endpoint of the interval.
    tol : float
        Tolerance for the stopping criterion.
    max_iter : int
        Maximum number of iterations.
    Returns
    -------
    x : float
        Approximate solution.
    """
    # creat a list to store the values of y
    y = []
    fa = f(a)
    fb = f(b)
    if fa*fb > 0:
        raise ValueError("f(a) and f(b) must have opposite signs.")
    for i in range(max_iter):
        x = (a*fb - b*fa)/(fb - fa)
        fx = f(x)
        y.append(x)
        if abs(fx) < tol:
            return y
        if fx*fa < 0:
            b = x
            fb = fx
        else:
            a = x
            fa = fx
    raise RuntimeError("Maximum number of iterations exceeded.")

### checking the number of iterations
def Table(x,y,func,e):
    
    print("Iteration     Bisection     Regula-Falsi")
    z = max(len(x),len(y))
    w = min(len(x),len(y))
    for i in range(z):
     if i < w:
      print("     ",i+1,"      ",x[i],"      ",y[i]) 
     else:
        if z==len(x):
            print("     ",i+1,"      ",x[i])
        else:
            print("      ",i+1,"         ",y[i])

##### Newton-Raphson method of root finding

def newton_raphson(f, df, x0, tol=1e-6, maxiter=100):
    """Newton-Raphson method of root finding.

    Args:
        f (_function_): Given function that is continous and differentiable.
        df (function): derivative of f.
        x0 (float): initial guess.
        tol (float): tolerance to check the convergence , Defaults to 1e-6.
        maxiter (int): maximum no of ineration , Defaults to 100.

    Raises:
        ValueError: "if the derivative of f is zero at x0 or if the function is not continous
        or does not satisfy convergence condition"

    Returns:
        float : root of the given function f.
    """
    # intial guess
    x = x0
    # iterate until convergence 
    for i in range(maxiter):
        xnew = x - f(x)/df(x)
        if abs(xnew - x) < tol:
            return xnew
        x = xnew
    # if no convergence bcz of bad guess or bad function, return None
    raise ValueError("Failed to converge after %d iterations" % maxiter)


#### laguerre and synthetic division method of root finding  

#To create fucntion for laguerre method
def Makefunc(c,j):
    n = len(c)
    """Find the value of the function at a given point.

    Args:
        c (array): co-efficients of the polynomial.
        j (float): point at which the function is to be evaluated.

    Returns:
        float: Value of the function at the given point.
    """
    y = 0.0
    for i in range(n):
        y = y + c[i]*(j**(n-i-1))
    return y
            
#For Differentiation of the function
 
def Diff(c):
    """co-efficients of the derivative of the polynomial.

    Args:
        c (Array): list of co-efficients of the polynomial.

    Returns:
        array: co-efficients of the derivative of the polynomial.
    """
    n = len(c)
    # create an empty list to store the co-efficients of the derivative
    m = []
    for i in range(n-1):
        m.append(c[i]*(n-1-i))
    return m
        
      
#To deflation of the function
 
def deflate(c,b):
    """Deflation of the polynomial.

    Args:
        c (array): co-efficients of the polynomial.
        b (float): root of the polynomial.

    Returns:
        array: quotient of the polynomial.
    """
    # create an empty list to store the co-efficients of the quotient
    m = []
    m.append(c[0])
    for i in range(1,len(c)-1):
        m.append(c[i]+b*m[i-1])
    return m  
    
#Laguerre method of root finding  

def Laguerre(c,b,e):
    """laguerre method of root finding.

    Args:
        b (float): initial guess.
        c (array): co-efficients of the polynomial.
        v(float): tolerance.
    Returns:
       float : root of the polynomial.
    """
    n = len(c)
    l = r = 0
    j = 0
    G,H = 0,0
    
    if Makefunc(c,b) == 0:
        return b
    # iterate until the root is found  
    while abs(b-l)>10**(-e):
        j = j+1
        l = b
        G = Makefunc(Diff(c),b)/Makefunc(c,b)
        H = G**2-Makefunc(Diff(Diff(c)),b)/Makefunc(c,b)
        r = math.sqrt((n-1)*(n*H-G**2))
      
        if abs(G-r) > abs(G+r):
            b = b-(n/(G-r))
        else:
            b = b-(n/(G+r))

    # print("Root:",b,", Iteration: ",j,"\n")
    return b
            
# solving the polynomial using laguerre method
def Solve(c,b,e):
    g=[]
    while len(c)>2:
        b = Laguerre(c,b,e)
        g.append(b)
        c = deflate(c,b)  
    g.append(-c[1]/c[0]) 
    print("Roots are:",g)    
    
    
#### Least square fit for any polynomial of degree n

def least_square_method(x, y, n):
    # n is the degree of the polynomial
    # x is the list of x values
    # y is the list of y values
    # x and y must have the same length
    # x and y must be lists
    # n must be an integer
    # n must be less than the length of x and y
    # n must be greater than 0
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n <= 0:
        raise ValueError("n must be greater than 0")
    if n >= len(x):
        raise ValueError("n must be less than the length of x and y")
    if not isinstance(x, list):
        raise TypeError("x must be a list")
    if not isinstance(y, list):
        raise TypeError("y must be a list")
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    # initializing the matrix
    matrix = []
    for i in range(n+1):
        matrix.append([])
        for j in range(n+1):
            matrix[i].append(0)
    # initializing the vector
    vector = []
    for i in range(n+1):
        vector.append(0)
    # calculating the matrix and vector
    for i in range(n+1):
        for j in range(n+1):
            for k in range(len(x)):
                matrix[i][j] += x[k]**(i+j)
        for k in range(len(x)):
            vector[i] += y[k]*x[k]**i
    # solving the matrix equation
    return a3.Gauss_Jordan(matrix, vector)

    