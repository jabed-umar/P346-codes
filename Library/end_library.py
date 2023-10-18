import math 
import matplotlib.pyplot as plt

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

### --------------- Random Number Generator by LCG method
# defing the class
class Random:
    def __init__(self, seed, range=[0,1]):
        self.seed = seed
        self.scale = lambda x: range[0] + x*(range[1]-range[0])
 # defining the LCG function
    def LCG(self, a = 1103515245, c = 12345, m = 2**32):
        """_summary_: Linear Congruential Generator

        Args:
            a (int, optional): Defaults to 1103515245.
            c (int, optional): Defaults to 12345.
            m (int, optional): Defaults to 2**32.

        Returns:
           list : list of random numbers with prefered range
        """ 
        self.seed = (a * self.seed + c) % m
        return self.scale(self.seed/m) 
    
## ------------ Random walk 2d 
# function to generate random numbers
def Random_bw(r0: float, n: int, k: int):
  a = 1103515245
  c = 12345
  m = 32768
  l = []

  for i in range (0, n):
    r0 = float(((a*r0 + c) % m)/m)
    p0 = 2*k*r0 - k
    l.append(p0)
  return l

#### --------- Defining a function to simulate random walk with given no. of steps

def random_walk(a: float, b: float, n: int):
  """_summary_: This function simulates a random walk with given no. of steps

  Args:
      a (float): seeds for random number generator in x direction
      b (float): seeds for random number generator in y direction
      n (int): no of steps
  """
  x0 = Random_bw(a, n, 1)
  y0 = Random_bw(b, n, 1)
  x1 = 0
  y1 = 0
  x = []
  y = []
  sum = 0
 
  for i in range(n):
    x1 += x0[i]
    x.append(x1)
    y1 += y0[i]
    y.append(y1)
    r = (x0[i]**2 + y0[i]**2)
    sum += r
  plt.plot(x, y)
  plt.title("Random walk with entered n steps")
  plt.show()
  rms = math.sqrt(sum/n)
  disp = math.sqrt(x0[n-1]**2 + y0[n-1]**2)
  print("rms distance for the walk = ", rms)
  print("net displacement for the walk = ", disp)


###------------------------Gauss Jordan method to solve linear equation 
## function to convert array to list 
def array_to_list(A,b):
    c = []
    for i in range(0,len(A)):
        c.append(b[0][i])
    return c
## --------------- Gauss_Jordan function
def Gauss_Jordan(A,b): 
    """_summary_ : Gauss Jordan method to solve linear equation Ax=b

    Args:
        A (2-array): Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-d array): Vector of order n made up of constants

    Returns:
        1d array: Solution vector x
    """
    # finding the number of rows/columns of the matrix A
    n = len(A)
    # finding the maximum element in each column
    for i in range(n):
        max = A[i][i]
        max_row = i
        for j in range(i+1,n):
            if abs(A[j][i]) > abs(max):
                max = A[j][i]
                max_row = j
        # swapping the row with the maximum element with the current row
        for k in range(n):
            A[i][k], A[max_row][k] = A[max_row][k], A[i][k]
        b[i], b[max_row] = b[max_row], b[i]
        # dividing the current row by the maximum element
        for l in range(n):
            A[i][l] = A[i][l]/max
        b[i] = b[i]/max
        # subtracting the current row from the other rows
        for m in range(n):
            if m != i:
                c = A[m][i]
                for o in range(n):
                    A[m][o] = A[m][o] - c*A[i][o]
                b[m] = b[m] - c*b[i]
    # printing the Row reduced (A) matrix 
    # print("The reduced A matrix is ", A)
    # printing the solution
    for p in range(n):
        b[p] = b[p]
    return b

### ------------------- solve a system of linear equation using doolittle decompose method 
def LU_decompose(A,b):
    """_summary_ : LU decomposition method to solve linear equation Ax=b

    Args:
        A (2-array): Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-d array): Vector of order n made up of constants

    Returns:
        1d array: Solution vector x
    """
    #finding the number of rows/columns
    n  = len(A)  
    #convert the matrix to upper and lower triangular matrix
    for j in range(n):
        for i in range(n):
            if i <= j :
                    sum = 0
                    for k in range(i):
                        sum += A[i][k]*A[k][j]
                    A[i][j] = A[i][j] - sum
            else  :
                    sum = 0
                    for k in range(j):
                        sum += A[i][k]*A[k][j]
                    A[i][j] = (A[i][j] - sum)/A[j][j]       
#forward substitution
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += A[i][j]*b[j]
        b[i] = b [i] - sum       
#backward substitution
    for i in range(n-1,-1,-1):
        sum = 0 
        for j in range(i+1,n):
            sum += A[i][j]*b[j]
        b[i] = (b[i] - sum)/(A[i][i])
    return b 


### ---------------- solve a system of linear equation using cholesky decompose method
# Check the given matrix is symmetric or not 
def Symmetric(mat):
    """_summary_ : Check if a matrix is symmetric

    Args:
        mat (2-d array): sqaure matrix of order n (n>=2)

    Returns:
        bool : True if the matrix is symmetric, else False
    """
    N = len(mat)
    for i in range(N):
        for j in range(N):
            if (mat[i][j] != mat[j][i]):
                return False
    return True
## Chelosky Decomposition function
def Decompose(A):
  # finding the number of rows/columns
 n = len(A)
 # decmpose into upper trangular matrix
 sum = 0
 for i in range(n):
  for j in range(i,n): 
   if j == i:
    for k in range(0,i):
     sum = sum + A[k][i]**2
    A[j][j]= round((A[j][j]-sum)**(0.5),4)   # rounding off to 4 decimal places
   else:
    A[j][i] = 0
    for k in range(0,i): 
     sum=sum+A[k][i]*A[k][j]
    A[i][j]= round((A[i][j]-sum)/A[i][i],4)      
   sum=0
 return A 
#Function for solving linear equations using Chelosky Decomposition     
def cholesky(A,b):
  # checking if the matrix is symmetric
  if Symmetric(A) == False:
    print("Matrix is not symmetric")
    return 0
  #Forward Substitution
  sum = 0
  A = Decompose(A)
  n = len(A) 
  for i in range(n):
    for k in range(i):
      sum=sum + A[k][i]*b[k][0]  # Taking transpose 
    b[i][0] = round((b[i][0]-sum)/A[i][i],4)    # rounding off to 4 decimal places
    sum=0
  #Backward Substitution 
  for i in range(n-1,-1,-1):
    for k in range(i+1,n):
        sum = sum + b[k][0]*A[i][k]  
    b[i][0] = round((b[i][0]-sum)/A[i][i],4)   # rounding off to 4 decimal places
    sum = 0
  return b

####--------- solve a system of linear equation using gauss_seidel method without rearranging the matrix
def Gauss_seidel(A,b,e):
    """_summary_ : Solve a system of linear equation using Gauss Seidal method

    Args:
        A (2-array):  Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-array): Vector of order n made up of constants
        e (precision): convergence criteria

    Returns:
       1_d array: Solution vector x
    """
    n = len(A)
    x = [[0] for y in range(n)]       # x stores values after new iteration
    m = 300
    y,sum = 0,0   
    for v in range(m):
     y = 0
     for i in range(n):
        sum = 0
        for j in range(n): 
            if j!= i:
             sum += A[i][j]*x[j][0] 
        c=(b[i][0]-sum)/A[i][i]
        if abs(c-x[i][0]) < (10**(-e)):    #Precision condition
            y += 1 
        x[i][0] = c 
     if y == n:   # If all elements of x follow precision condition
         break  
    print("Number of iterations is :", v+1,"\nThe solution matrix x is:\n") 
    print(x)

###--------------------solve a system of linear equation using doolittle decompose method 
def LU_decompose(A,b):
    """_summary_ : LU decomposition method to solve linear equation Ax=b

    Args:
        A (2-array): Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-d array): Vector of order n made up of constants

    Returns:
        1d array: Solution vector x
    """
    #finding the number of rows/columns
    n  = len(A)  
    #convert the matrix to upper and lower triangular matrix
    for j in range(n):
        for i in range(n):
            if i <= j :
                    sum = 0
                    for k in range(i):
                        sum += A[i][k]*A[k][j]
                    A[i][j] = A[i][j] - sum
            else  :
                    sum = 0
                    for k in range(j):
                        sum += A[i][k]*A[k][j]
                    A[i][j] = (A[i][j] - sum)/A[j][j]       
#forward substitution
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += A[i][j]*b[j]
        b[i] = b [i] - sum       
#backward substitution
    for i in range(n-1,-1,-1):
        sum = 0 
        for j in range(i+1,n):
            sum += A[i][j]*b[j]
        b[i] = (b[i] - sum)/(A[i][i])
    return b 


###------------solve a system of linear equation using Guass_Seidel method with rearranging the matrix
#Pivoting function for producing diagonally dominant matrix
def pivot(A):
    """_summary_ : Pivoting function for producing diagonally dominant matrix

    Args:
        A (2-d array): sqaure matrix of order n (n>=2)

    Returns:
        2-d array : Diagonally dominant matrix
    """
    c = 0
    t = 0
    for i in range(len(A)):   #Row pivot
     t = i                 # t stores largest element of a column
     c = abs(A[i][i])      # taking absolute value of diagonal elements 
     for j in range(len(A[0])):
         if abs(A[i][j]) > c:   # checking if the element is greater than the diagonal element
             c = abs(A[i][j])
             t = j              
     if t > i:
         for k in range(len(A)):
             A[k][i],A[k][t]= A[k][t],A[k][i]   
     elif t < i:
         print("Matrix is not diagonally dominant \n")
         return 0  
    return A


#----------------------------Jacobi and Guass-Seidal method to solve the linear equation 
## Jacobi function     
def Jacobi(A,b,e):
    """_summary_ : Solve a system of linear equation using Jacobi method

    Args:
        A (2-array):  Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-array): Vector of order n made up of constants
        e (precision): convergence criteria
    Returns:
         1_d array: Solution vector x
    """
    A = pivot(A)
    if A == 0:  #From pivot function
      print("Jacobi not possible")
      return 0 
    n = len(A)
    C = [[1] for y in range(n)]       # C stores values after new iteration
    D = [[0] for y in range(n)]        # D stores the values after last iteration
    m = 50000                        # m stores maximum number of iterations
    sum = 0
    y = 1
    for k in range(m):
     for i in range(n):
        for j in range(n):
            if j != i:
             sum = sum+A[i][j]*C[j][0]
            if abs(D[j][0]-C[j][0]) > (10**(-e)):y=1  #Checking for precision        
        if y==1:    
         D[i][0]=(b[i][0]-sum)/A[i][i]
        else:
            break
        sum=0  
     y = 0    
     C,D = D,C 
    print("Number of iterations is:",k+1,"\nThe solution matrix x is:\n")
    print(C)   
    
    
##----------------Guass Seidel function to solve linear equations
def Gauss_Seidel(A,b,e):
    """_summary_ : Solve a system of linear equation using Gauss Seidal method

    Args:
        A (2-array):  Square matrix of order n (n>=2) made up of coefficinets of the variables
        b (1-array): Vector of order n made up of constants
        e (precision): convergence criteria

    Returns:
       1_d array: Solution vector x
    """
    n = len(A)
    A = pivot(A)
    if A == 0: # From pivot function checking the diagoanl dominance
      print("Guass-Seidal not possible")
      return 0
    x = [[0] for y in range(n)]       # x stores values after new iteration
    m = 300
    y,sum = 0,0   
    for v in range(m):
     y = 0
     for i in range(n):
        sum = 0
        for j in range(n): 
            if j!= i:
             sum += A[i][j]*x[j][0] 
        c=(b[i][0]-sum)/A[i][i]
        if abs(c-x[i][0]) < (10**(-e)):    #Precision condition
            y += 1 
        x[i][0] = c 
     if y == n:   # If all elements of x follow precision condition
         break  
    print("Number of iterations is:", v+1,"\nThe solution vector x is:\n") 
    print(x)
    
##-------------------- Transpose of a matrix
def transpose(A):
    """This function calculates the transpose of a matrix
    Args:
        A (list): This is the matrix
    """
    n = len(A)
    m = len(A[0])
    B = [[0 for i in range(n)] for j in range(m)]
    for i in range(n):
        for j in range(m):
            B[j][i] = A[i][j]
    return B
    
### Importing the libraries
import math
import Library.assign3 as a3


#----------------------------------------Bracket function to find the interval
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
      
      
#### -----------------------------------------bisection method to find the root of a function

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

####----------------------------------------------------- Regula falsi method of root finding

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

##### ----------------------------------------Newton-Raphson method of root finding

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


####--------------------------------laguerre and synthetic division method of root finding  

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
    print("Roots of the polynomial is :",g)    
    
    
####------------------------ Least square fit for any polynomial of degree n

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
    return a3.Gauss_Jordan(matrix, vector)

###---------------------------------------- Calculate pearson correlation coefficient
def pearson(x,y):
    """This function calculates the pearson correlation coefficient
    Args:
        x (list): list of x values
        y (list): list of y values
    """
    n = len(x)
    sum_x = 0
    sum_y = 0
    sum_xy = 0
    sum_x2 = 0
    sum_y2 = 0
    for i in range(n):
        sum_x += x[i]
        sum_y += y[i]
        sum_xy += x[i]*y[i]
        sum_x2 += x[i]**2
        sum_y2 += y[i]**2
    # calculate the pearson correlation coefficient
    r = (n*sum_xy - sum_x*sum_y)/((n*sum_x2 - sum_x**2)**0.5*(n*sum_y2 - sum_y**2)**0.5)
    print("The pearson correlation coefficient is ", r)
    return r

###------------------------------------- integration by mid point method 
def midpoint(f, a, b, n):
    h = (b-a)/n
    s = 0
    for i in range(n):
        s += f(a + h*(i+0.5))
    # print ("For n =", n)
    return s*h


###--------------------------------------integration by trapezoidal method
def trapezoidal(f, a, b, n):
    h = (b-a)/n
    s = 0.5*(f(a) + f(b))
    for i in range(1, n):
        s += f(a + h*i)
    print ("For n =", n)
    return s*h


###---------------------------------Integration by Simpson Method (1)
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

###---------------------- integration by simpson's method (2)
def Simpson(funx,a,b,N):
    k=1
    sum=0
    h=(b-a)/N
    
    for i in range(N):
        if i==0 and i==N-1:
            k=1
        elif i%2==0:
            k=2
        else:
            k=4
        sum= sum + k*funx(a+i*h)

    sum= round(sum*h/3, 8)
    # print("sum for N (",N,") =",sum)
    # print ("For n =",N)
    return sum

### integration by Monte Carlo method
import random
def monte_carlo(f, a, b, n):
    s = 0
    for i in range(n):
        x = random.uniform(a, b)
        s += f(x)
    print ("For n =",n)
    return s*(b-a)/n

## define the function for the linear mass density function
def rho(x):
    return x**2

# ------------------------ Finding the value of N for which the error is less than 0.0001
# first function for mid point method 
def max_N(a,b,f,ddf,e):
    return math.ceil(math.sqrt((b-a)**3*abs(ddf)/(24*e)))
# 2nd function for trapezoidal method
def max_n(a,b,fun,ddf,e):
    return math.ceil(math.sqrt((b-a)**3*abs(ddf)/(12*e)))
# 3rd function for simpson's method
def max(a,b,fun,ddddf,e):
    return math.ceil(((b-a)**3*abs(ddddf)/(180*e))**(1/4))


# solving a differential equation 
# ----------------------------forwad euler method
def Euler_f(f, x0, y0, h, n):
    """Solve the differential equation by euler forward method y'=f(x,y) with initial condition y(x0) = y0

    Args:
        f (func): the given function f(x,y)
        x0 (float): the initial value of x
        y0 (float): the initial value of y
        h (float): the step size
        n (int): the number of steps

    Returns:
        x (float), y(float): the final value of x and y
    """
    # initialize x and y
    x = x0
    y = y0
    # create a list to store the values of x and y
    a = []
    b = []
    # main loop
    for i in range(n):
        y = y + h * f(x, y)
        a.append(y)
        x = x + h
        b.append(x)
        plt.plot(b, a)
    return x,y

#------------------------ corrector method
def predictor_corrector(f, x0, y0, h, n):
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
        k2 = h * f(x + h, y + k1)
        # calculate y
        y = y + (k1 + k2) / 2
        a.append(y)
        x = x + h
        b.append(x)
        plt.plot(b, a) 
    return x,y
#--------------------------------4th order Runge-Kutta method for solving differential eqn
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
    
### ------------------- Solving the differential equation using the Runge-Kutta method
# coupled differential equation
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


## ------------------------------ Radioactive decay
# Radioactive decay of a sample of N atoms of a radioactive isotope
import numpy as np
import matplotlib.pyplot as plt
# Half life of the isotope
t_half = 10
# decay constant
lam = np.log(2)/t_half
# time step
dt = 1
# decay probability
p = lam
# survival probability
q = 1 - p
# initial number of atoms
N = 1000
# define a function to calculate the number of atoms at time t
def decay(n):
    population = []
    for t in range(100):
        r = np.random.random(n)
        # number of atoms that survive
        survived = np.sum(r < q)
        population.append(survived)
        n = survived
    return population

# plot the population of atoms at different times
# plt.plot(range(100), decay(N))
mean_decay = np.mean(np.array([decay(N) for i in range(100)]), axis=0)
std_decay = np.std(np.array([decay(N) for i in range(100)]), axis=0)
T = range(100)
#plot the mean decay with time
#plt.plot(T, mean_decay)
exact = N*np.exp(-lam*T)
# plot the exact theoretical decay with time
#plt.plot(T, exact)

## ------------------------------ Radioactive decay with daughter nuclei
def decaya2b(na,thalf,dt,t):
    l1 = 0.693/thalf
    n = int(t/dt)
    import random
    import matplotlib.pyplot as plt
    import numpy as np
    nb = 0
    bx = []
    by = []
    ax = []
    ay = []
    for j in range(0,n):
        for i in range(0,na):
            if random.random()<= l1*dt :
                na = na - 1
                ay.append(na)
                ax.append((j+1)*dt)
                nb = nb + 1
                by.append(nb)

    xpoints = np.array(ax)
    ypoints = np.array(by)
    #plt.plot(xpoints, ypoints, 'o')
    xpoints = np.array(ax)
    ypoints = np.array(ay)
    #plt.plot(xpoints, ypoints)
    #plt.xlabel('X-axis')
    #plt.ylabel('Y-axis')
    #plt.title("Decay Plot")
    #plt.show()
    
    