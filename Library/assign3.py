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
#Function to make tarnspose of matrix
def transpose(A):
    B=[]
    for i in range(len(A)):
        B.append([])
    for j in range(len(A)):
        B[i].append(A[j][i])
    return B


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