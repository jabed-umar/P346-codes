###1(a). Gauss Jordan method to solve linear equation 
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

###1(b). solve a system of linear equation using doolittle decompose method 
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


### 2(a). solve a system of linear equation using cholesky decompose method
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

#### 2(b). solve a system of linear equation using gauss_seidel method without rearranging the matrix
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
    print("Number of iterations:", v+1,"\nSolution matrix X:\n") 
    print(x)

### 3(a). solve a system of linear equation using doolittle decompose method 
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


### 3(b). solve a system of linear equation using Guass_Seidel method with rearranging the matrix
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
    
