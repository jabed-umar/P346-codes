
### ------------------- Solving the differential equation using the Runge-Kutta method

def rk_4(d2ydx2, dydx, x0, y0, z0, xf, h):
    """Runge-Kutta 4th order method for solving 2nd order ODEs.

    Args:
        d2ydx2 (list): double derivative of y with respect to x.
        dydx (list): single derivative of y with respect to x.
        x0 (float): initial value of x.
        y0 (float): initial value of y.
        z0 (float): initial value of z.
        xf (float): final value of x.
        h (float): strip size.

    Returns:
        list: solution of the ODE.
    """
    x = [x0]
    y = [y0]
    # dy/dx = z 
    df = [z0]    
    # no. of steps to take 
    n = int((xf-x0)/h)     
    for i in range(n):
        x.append(x[i] + h)
        k1 = h * dydx(x[i], y[i], df[i])
        l1 = h * d2ydx2(x[i], y[i], df[i])
        k2 = h * dydx(x[i] + h/2, y[i] + k1/2, df[i] + l1/2)
        l2 = h * d2ydx2(x[i] + h/2, y[i] + k1/2, df[i] + l1/2)
        k3 = h * dydx(x[i] + h/2, y[i] + k2/2, df[i] + l2/2)
        l3 = h * d2ydx2(x[i] + h/2, y[i] + k2/2, df[i] + l2/2)
        k4 = h * dydx(x[i] + h, y[i] + k3, df[i] + l3)
        l4 = h * d2ydx2(x[i] + h, y[i] + k3, df[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        df.append(df[i] + (l1 + 2*l2 + 2*l3 + l4)/6)

    return x, y


### defining the matrix multiplication function
def matrix_product(A, B):
    """Matrix product of two matrices.

    Args:
        A (list): matrix A.
        B (list): matrix B.

    Returns:
        list: matrix product of A and B.
    """
    # no. of rows in A
    m = len(A)
    # no. of columns in A
    n = len(A[0])
    # no. of rows in B
    p = len(B)
    # no. of columns in B
    q = len(B[0])
    # check if matrix multiplication is possible
    if n != p:
        raise Exception ("Matrix multiplication is not possible")
    # initialize the product matrix
    C = [[0 for i in range(q)] for j in range(m)]
    for i in range(m):
        for j in range(q):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

#### Function to find dot product of two matrices
def dot(A,B):
    """Find dot product of two matrices.

    Args:
        A (array): matrix A.
        B (array): matrix B.

    Returns:
        float: dot product of A and B.
    """
    # initialize the sum
    sum = 0
    for i in range(len(A[0])):
        sum += A[0][i]*B[0][i]
    return sum

### ----------------------------- Function for returning dominant eigenvalue and corresponding eigenvector
def powr_iter(A,x0,e):
    """Find dominant eigenvalue and corresponding eigenvector.

    Args:
        A (array): a matrix.
        x0 (list): guess vector.
        e (float): tolerance.

    Returns:
        float & array: eigenvalues and eigenvectors of matrix A
    """
    sum = 2
    y = x0
    z = matrix_product(A,y)
    k0 = dot(z,x0)/dot(y,x0)
    y = z
    z = matrix_product(A,y)
    k1 = dot(z,x0)/dot(y,x0)
    
    while abs(k1-k0)> e:
      sum += 1
      k0 = k1
      y = z
      z = matrix_product(A,y)  
      k1 = dot(z,x0)/dot(y,x0)
    
    ####  normalized eigenvector
    sum=0 
    for i in range(len(y)):
        sum += y[i][0]**2
    sum = sum**0.5
    for i in range(len(y)):
        y[i][0] /= sum
    return y, k1

### --------------------------------- Solving the 1-d heat equation 
def heat_equation(temp0:callable, L: float, T:float, nL:int, nT:int, t_upto:int = None):
    """Solves the heat equation 

    Args:
        temp0 (callable): initial temperature distribution.
        L (float): length of the rod.
        T (float): time period.
        nL (int): no. of strips in the length.
        nT (int): no. of strips in the time.
        t_upto (int): time upto which the solution is to be plotted. default is nT.

    Returns:
        2D list: solution of the heat equation.
    """
    if t_upto is None: t_upto = nT
    
    ht = T/nT
    hx = L/nL
    alpha = ht/(hx**2)
    # print(alpha)
    
    A = [[0 for i in range(nL)] for j in range(t_upto)]
    for i in range(nL):
        A[0][i] = temp0(i, nL)
    
    for t in range(1, t_upto):
        for x in range(nL):
            if x == 0:
                A[t][x] = A[t-1][x]*(1-2*alpha) + A[t-1][x+1]*alpha
            elif x == nL-1:
                A[t][x] = A[t-1][x-1]*alpha + A[t-1][x]*(1-2*alpha)
            else:
                A[t][x] = A[t-1][x-1]*alpha + A[t-1][x]*(1-2*alpha) + A[t-1][x+1]*alpha
    return A

### ------------------------------ Shooting method
# RK for shooting method
def rk_shoot(d2ydx2, dydx, x0, y0, z0, xf, h):
    """shooting method using RK4.

    Args:
        d2ydx2 (array): derivative of the second order.
        dydx (array): derivative of the first order.
        x0 (float): initial value of x.
        y0 (float): initial value of y.
        z0 (float): initial value of z.
        xf (float): final value of x.
        h (float): strip size.

    Returns:
        array: solution of the differential equation.
    """
    x = [x0]
    y = [y0]
    z = [z0]
    # no of strips
    N = int((xf-x0)/h)
    for i in range(N):
        k1 = h * dydx(x[i], y[i], z[i])
        l1 = h * d2ydx2(x[i], y[i], z[i])

        k2 = h * dydx(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)
        l2 = h * d2ydx2(x[i] + h/2, y[i] + k1/2, z[i] + l1/2)

        k3 = h * dydx(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)
        l3 = h * d2ydx2(x[i] + h/2, y[i] + k2/2, z[i] + l2/2)

        k4 = h * dydx(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * d2ydx2(x[i] + h, y[i] + k3, z[i] + l3)

        x.append(x[i] + h)
        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (l1 + 2*l2 + 2*l3 + l4)/6)
    return x, y, z

# Lagrange interpolation for the intrpolation part 
def lag_inter(zeta_h, zeta_l, yh, yl, y):
    zeta = zeta_l + (zeta_h - zeta_l) * (y - yl)/(yh - yl)
    return zeta

# Shooting method for solving the 2nd order ODE
def shoot(d2ydx2, dydx, x0, y0, xf, yf, z1, z2, h, tol=1e-6):  
    x, y, z = rk_shoot(d2ydx2, dydx, x0, y0, z1, xf, h)
    yn = y[-1]
    if abs(yn - yf) > tol:
        if yn < yf:
            zeta_l = z1
            yl = yn
            x, y, z = rk_shoot(d2ydx2, dydx, x0, y0, z2, xf, h)
            yn = y[-1]
            if yn > yf:
                zeta_h = z2
                yh = yn
                zeta = lag_inter(zeta_h, zeta_l, yh, yl, yf)
                x, y, z = rk_shoot(
                    d2ydx2, dydx, x0, y0, zeta, xf, h)
                return x, y
            else:
                print("Invalid bracketing.")
        elif yn > yf:
            zeta_h = z1
            yh = yn
            x, y, z = rk_shoot(d2ydx2, dydx, x0, y0, z2, xf, h)
            yn = y[-1]
            if yn < yf:
                zeta_l = z2
                yl = yn
                zeta = lag_inter(zeta_h, zeta_l, yh, yl, yf)
                x, y, z = rk_shoot(
                    d2ydx2, dydx, x0, y0, zeta, xf, h)
                return x, y
            else:
                print("Invalid bracketing.")
    else:
        return x, y
