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
## Calculate the area of ellipse using Monte Carlo method
## Area of ellipse by random number 
def area(n):
    x = 2*My_random(10,n,0) # Change the range x into [0,2]
    y = My_random(10,n,0)
    inside = 0 
    # check if the point is inside the ellipse or not
    for i in range(n):
        if x[i]**2+ 4*y[i]**2 <= 4:
            inside += 1
    # calculate the volume of the ellipse
    vol_e = 8*inside/n
    return vol_e


###2. Newton's method for finding square root
def newton_raphson(f, df, x0, tol=1e-4, maxiter=100):
    """Newton-Raphson method of root finding.

    Args:
        f (_function_): Given function that is continous and differentiable.
        df (function): derivative of f.
        x0 (float): initial guess.
        tol (float): tolerance to check the convergence , Defaults to 1e-4.
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


##3. solve a system of linear equation using gauss_seidel method without rearranging the matrix
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
    

###4. Least square method with error = 0
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
    print("The slope of the linear fit is ", slope)
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
    return r