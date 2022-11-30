
### 3.------------------------ Newton Raphson method
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

### 6.----------------------------- Function for returning eigenvalues and corresponding eigenvectors
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

### Function for returning dominant eigenvalue and corresponding eigenvector
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
    sum = 0 
    for i in range(len(y)):
        sum += y[i][0]**2
    sum = sum**0.5
    for i in range(len(y)):
        y[i][0] /= sum
    return y, k1

### 4. ---------------------------------Integration by Simpson Method (1)
def simpson(f, x, y, n,flag):
        # calculating n 
        #flag is used to gice output depending on requirement 
        h = abs(y - x) / n
        sn = f(x) + f(y)
        # applying the algorithm
        for i in range(1, n):
            if i % 2 == 0:
                x = x + h
                sn = sn + 2 * f(x)
            else:
                x = x + h
                sn = sn + 4 * f(x)
        t = (h / 3) * sn
        if flag ==True :
            print("Integration of the given function for n = ", n ,"is :",round(t,9))
        else:
            return round(t,9)
        
#### 7. ------------------------ Least square fit for any polynomial of degree n
import Library.assign3 as assign3
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
    # solving the matrix equation to find the coefficients of the polynomial
    # calling the gaussian elimination method to calculate the coefficients
    return assign3.Gauss_Jordan(matrix, vector)     


#5.--------------------------------4th order Runge-Kutta method for solving differential eqn
# uncoupled differential equation
def range_kutta_fourth(f, x0, y0, h, n):
    # initialize x and y
    x = x0
    y = y0
    # create a list to store the values of x and y
    a = []
    b = []
    # calculate k1
    for i in range(n):
        k1 = h * f(x, y)
        # calculate k2
        k2 = h * f(x + h / 2, y + k1 / 2)
        # calculate k3
        k3 = h * f(x + h / 2, y + k2 / 2)
        # calculate k4
        k4 = h * f(x + h, y + k3)
        # calculate y
        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        a.append(y)
        x = x + h
        b.append(x)
        plt.plot(b, a)
    return x,y
## 2. --------------------------- Function for the equilibria of the system
### ----------- Random Number Generator by LCG method
import matplotlib.pyplot as plt

class Random:
    def __init__(self, seed: float = 0.1, range: list = [0,1]):
        """Random number generator.
        Args:
            seed (float, optional): Initial Seed. Defaults to 0.1.
        """
        self.seed = seed
        self.scale = lambda x: range[0] + x*(range[1]-range[0])

    def rand(self, c = 3.5):
        """A function to generate a random number by using the formula:
            X[n+1] = c * X[n] * (1- X[n]).
            
        Args:
            c (float, optional): Defaults to 3.5.
            
        Returns:
            float: a random number between 0 and 1.
        """
        self.seed = c*self.seed*(1-self.seed)
        return self.seed

    def LCG(self, a = 1103515245, c = 12345, m = 32768):  # LCG
        self.seed = (a*self.seed + c) % m
        return self.scale(self.seed / m)
    
    
# 1.  ---------------------------------Lu decomposition of given matrix
def LU(A):
    n = len(A)
    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]
    for i in range(n):
        L[i][i] = 1.0
        for j in range(i, n):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = A[i][j] - s1
        for j in range(i, n):
            if i == j:
                L[i][i] = 1.0
            else:
                s2 = sum(U[k][i] * L[j][k] for k in range(i))
                L[j][i] = (A[j][i] - s2) / U[i][i]
    return L, U

# checking the inverse of the matrix 
# determinant of a given matrix 
def determinant(A):
    # LU decomposition
    L, U = LU(A)
    # determinant of A
    det = 1
    for i in range(len(A)):
        det *= L[i][i] * U[i][i]
    return det

# inverse of a given matrix
def inverse(A):
    # LU decomposition
    L, U = LU(A)
    # inverse of A
    inv = [[0.0] * len(A) for i in range(len(A))]
    for i in range(len(A)):
        inv[i][i] = 1.0 / L[i][i]
        for j in range(i):
            s = sum(U[k][j] * inv[i][k] for k in range(j))
            inv[i][j] = -s / L[j][j]
        for j in range(i + 1, len(A)):
            s = sum(U[k][j] * inv[i][k] for k in range(i))
            inv[i][j] = -s / L[i][i]
    return inv