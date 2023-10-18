### Random Number Generator by LCG method
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


### volume of a sphere by random number
def vol(n):
    """This function calculates the volume of a sphere using random number
    Args:
        n (integer): This is the number of random numbers in the list
    """
    x = My_random(10,n,0)
    y = My_random(5,n,0)
    z = My_random(3,n,0)
    inside = 0
    # check if the point is inside the sphere
    for i in range(n):
        if x[i]**2 + y[i]**2 + z[i]**2 <= 1:
            inside += 1
    # calculate the volume of the sphere
    vol_s = inside/n
    print(vol_s)
    return vol_s
    # print(len(x))
    
### Calculate the value of pi using Monte Carlo method
def pi_value(n):
    """This function calculates the value of pi using Monte Carlo method
    Args:
        n (integer): This is the number of random numbers in the list
    """
    x = My_random(10,n,0)
    y = My_random(5,n,0)
    inside = 0
    # check if the point is inside the circle
    for i in range(n):
        if x[i]**2 + y[i]**2 <= 1:
            inside += 1
    pi = 4*inside/n
    print(pi)
    return pi

### Newton's method for finding square root
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
 
 
### evaluate the derivative of a function at a point x by secant method
def derivative(f, x, h=1e-6):
    """This function evaluates the derivative of a function at a point x
    Args:
        f (function): The function whose derivative is to be evaluated
        x (float): The point at which the derivative is to be evaluated
        h (float, optional): The step size. Defaults to 1e-6.
    Returns:
        float: The derivative of the function at the point x
    """
    return (f(x+h) - f(x-h))/(2*h)


###solve a system of linear equation using Guass_Seidel method with rearranging the matrix
# Transpose of a matrix
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
#Jacobi and Guass-Seidal method to solve the linear equation 
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
## Guass Seidel function to solve linear equations
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


### Lagrange interpolation
def interpolation(x, y, x0):
    """Interpolate y(x) at x0 using Lagrange polynomial."""
    n = len(x)
    y0 = 0.0
    for i in range(n):
        p = 1.0
        for j in range(n):
            if i != j:
                p *= (x0 - x[j])/(x[i] - x[j])   # Lagrange polynomial
        y0 += p*y[i]
    return y0

### Least square method with error = 0
def least_sqr_method(x,y):
    """This function calculates the least square method
    Args:
        x (list): list of x values
        y (list): list of y values
    """
    n = len(x)
    sum_x = 0
    sum_y = 0
    sum_xy = 0
    sum_x2 = 0
    for i in range(n):
        sum_x += x[i]
        sum_y += y[i]
        sum_xy += x[i]*y[i]
        sum_x2 += x[i]**2
    # calculate the slope and intercept
    slope = (n*sum_xy - sum_x*sum_y)/(n*sum_x2 - sum_x**2)
    intercept = (sum_y - slope*sum_x)/n
    print("The slope is ", slope)
    print("The intercept is ", intercept)
    return slope, intercept

## Calculate pearson correlation coefficient
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


### Least sqr fit with error
import numpy as np  
import matplotlib.pyplot as plt  

def graph(formula, x_range):  
    x = np.array(x_range)  
    y = formula(x)  # <- note now we're calling the function 'formula' with x
    plt.plot(x, y)  
    plt.show()  
    
def Fit(A,B):
 n=3
 sum1=0
 sum2=0
 C=[[0 for x in range(n)] for y in range(n)]
 D=[[0] for y in range(n)]
 plt.plot(A,B)
 
 for i in range(3):
     for j in range(i+1):
             
         for k in range(len(A)):
           sum1=sum1+(A[k][0]**(i+j))
         
         C[i][j]=C[j][i]=sum1
         sum1=0
     for k in range(len(A)):
        sum2=sum2+((A[k][0]**(i))*B[k][0])  
     D[i][0]=sum2
     sum2=0
 
 
 for line in C:
   print ('  '.join(map(str, line))) 
 print()
 for line in D:
   print ('  '.join(map(str, line))) 
 print()
 return C,D

 
     

def my_formula(x): return 4.004+2.393*x-0.029*(x**2)       

graph(my_formula, range(20,120))
