import matplotlib.pyplot as plt
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
def My_random(seed,n,k):
    """This function is used to generate random number by LCG method
    Args:
        seed (integer): This is the seed for the random number generator
        n (integer): Number of the of the random number 
        k (Either 0 or anything): Just to change the range of random numbers
    Returns:
        Floats : list of random numbers
    """
    a = 1103515245 
    c = 12345 
    m = 32768
    x = seed
    rand= []                # This is the list of random numbers
    # 100 random numbers
    if k == 0:
     for i in range(n):    
    # LCG method for creating random number in the range [0,1]
        x = (a * x + c) % m 
        x = x/m      # Normalizing the random number
        rand.append(x)
    else:
     for i in range(n):    
    # LCG method for creating random number in the range [-1,1] 
        x = 2*((a * x + c) % m)/m -1     # Normalizing the random number
        rand.append(x)   
    return rand
   
## 1.____________________________________________Finding inverse by LU decompose method
# forward and backward substitution
def for_back_subs(mat_L, mat_U, vec_B):
    #initialize two vectors to store the solution
    x = [0] * len(mat_U[0])
    y = [0] * len(mat_L[0])
    
    #forward substitution    
    for i in range(len(mat_L)):
        sum = 0
        for j in range(i):
            sum += mat_L[i][j]*y[j]
        y[i] = vec_B[i] - sum

    #backward substitution (loop should run in reverse)
    for i in reversed(range(len(mat_U))):
        sum = 0
        for j in reversed(range(i, len(mat_U[0]))):
            sum += mat_U[i][j]*x[j]
        x[i] = (y[i] - sum)/mat_U[i][i]
    return x
# LU decompose method
def LU(A):
    """_summary_: LU decomposition of a matrix

    Args:
        mat_A (array): given matrix

    Returns:
        arrys: L and U matrices
    """
    n = len(A)
    #initialize empty matrix to store decomposed matrices
    L = [[0 for x in range(n)] for y in range(n)] 
    U = [[0 for x in range(n)] for y in range(n)]

    #loop through every columns of A 
    for i in range(len(A)):
        #Diagonal entries of L matrix are all 1
        L[i][i] = 1

        #Upper triangular matrix
        for k in range(i, len(A[0])):
            sum = 0
            for j in range(i): 
                sum += (L[i][j] * U[j][k])
  
            U[i][k] = A[i][k] - sum

        #Lower triangular matrix
        for k in range(i+1, len(A[0])):             
            sum = 0 
            for j in range(i): 
                sum += (L[k][j] * U[j][i]) 
            
            #The order of element is reversed from i,k to k,i
            L[k][i] = (A[k][i] - sum) / U[i][i] 
  
    return L,U
# inverse of a matrix cheking by the determinant
def check_inverse(C):
    """_summary_: Check if the inverse of a matrix exist by determinant checking

    Args:
        mat_U (array): given matrix

    Raises:
        Warning: if determinant is 0

    Returns:
        Bool: If the inverse exist return True, else return False
    """
    det = 1
    for i in range(len(C)):
        det = det*C[i][i]

    if det == 0:
        raise Warning("Inverse of the matrix does not exist")
        return False
    
    else:
        return True
# finding the inverse of a matrix
def matrix_inverse(mat_A):
    """_summary_: Inverse of a matrix

    Args:
        mat_A (array): given matrix 

    Returns:
        array: inverse of the given matrix
    """
    # check the no of rows and columns
    n = len(mat_A)
    # initialise the inverse matrix with 0
    inverse_matrix = [[0 for x in range(n)] for y in range(n)] 
    # calling LU decomposition function to decompose A into L and U
    L,U = LU(mat_A)
    #Check if the inverse exist
    if check_inverse(U) == True:
        for i in range(n):
            #initialize vector with ith value 1, rest 0
            B = [0 for x in range(n)]
            B[i] = 1
            column = for_back_subs(L, U, B)
            #Store the solution column wise
            for j in range(n):
                inverse_matrix[j][i] = round(column[j], 3)
        return inverse_matrix

    else:
        return None    
