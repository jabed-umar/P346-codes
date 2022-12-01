import numpy as np
import matplotlib.pyplot as plt
# ----------------------------------------performing SVD by using numpy
def SVD(A):
    """_summary_: This function performs SVD on the given matrix

    Args:
        A (array): The matrix on which SVD is to be performed

    Returns:
        arrys: It returns the U, S and Vt matrices
    """
    U, S, Vt = np.linalg.svd(A)    
    return U, S, Vt


# ------------------------------------------ Contruct the sigma matrix 
def sigma_matrix(A):
    """_summary_: This function constructs the sigma matrix of A(m*n) matrix

    Args:
        S (array): The singular values array

    Returns:
        array: It returns the sigma matrix
    """
    # calling the SVD function
    U, S, Vt = np.linalg.svd(A)
    
    # constructing the sigma matrix for square matrix
    if len(A) == len(A[0]):
        n = np.diag(S)
    
    # constructing the sigma matrix for rectangular matrix with m > n
    elif len(A) > len(A[0]):
        r = (len(U)-len(Vt))
        c = np.zeros((r,len(Vt)))
        n = np.vstack([np.diag(S), c])
        
    # constructing the sigma matrix for rectangular matrix with m < n
    else:
        r = (len(Vt)-len(U))
        c = np.zeros((len(U),r))
        n = np.hstack([np.diag(S), c])
    return n

## ------------------------------------------- Low rank approximation function
# Changing the U and Vt matrix to the for the image compression
def compress_matrix(n,B):
    U,D,Vt = SVD(B)
    k = len(np.diag(D))
    if n > k:
        raise ValueError("n must be less than the number of singular values")
    elif n <= 0:
        raise ValueError("n must be greater than 0")
    elif n == k:
        return U, np.diag(D), Vt
    else:
        U = np.delete(U, np.s_[n:], axis=1)
        D = np.delete(D, np.s_[n:])
        Vt = np.delete(Vt, np.s_[n:], axis=0)
        return U, np.diag(D), Vt


# ---------------------------------- function to find the transpose of a matrix 

def transpose(a):
    if len(a) != len(a[0]):
        raise Exception ("Matrix is not square")
    else:
        return [[a[j][i] for j in range(len(a))] for i in range(len(a[0]))]


#------------------------------------------- defining the matrix product function
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


## ------------------------------------------- variance of singular values function
def variance(X):
    """_summary_ : This function is used to find the variance of the singular values

    Args:
        X (array): The matrix of the image

    Returns:
        array: It returns the variance of the singular values
    """
    x,y,z = SVD(X)
    return np.round(y**2/np.sum(y**2),6)


## ------------------------------------------- image compression function for gray scale image
#--------- with numpy inbuilt matrix multiplication function
def compress_grey(rank,A):
    """_summary_ : This function is used to compress the image by using the rank of the image

    Args:
        rank (int): The rank of the image
        A (array): The matrix of the image
    """
    U,s,V = SVD(A)
    compress_a = U[:, :rank]@ np.diag(s[:rank])@V[:rank, :]
    compress_a = compress_a.astype(int)
    return compress_a
    # plt.imshow(compress_a, cmap = "gray")


# ------------------------------------------- image compression function for color 
# with numpy inbuilt matrix multiplication function
def compress_color_a(rank,R,G,B):
    """_summary_ : This function is used to compress the image by using the rank of the image

    Args:
        rank (int): The rank of the image
        R (array ): The matrix of the image because of red channel  
        B (array): The matrix of the image because of blue channel
        G (array): The matrix of the image because of green channel
    """
    # calling the SVD function
    U, S, Vt = SVD(R)  # for red channel
    u,s,Vt = SVD(G)   # for green channel
    u,S,Vt = SVD(B)   # for blue channel
    
    
    compress_red = U[:, :rank]@ np.diag(S[:rank])@Vt[:rank, :]   # @ is used for matrix multiplication
    compress_red = compress_red.astype(int)
    #print(compress_red)
    compress_green = u[:, :rank]@ np.diag(s[:rank])@Vt[:rank, :]  # @ is used for matrix multiplication
    compress_green = compress_green.astype(int)
    #print(compress_green)
    compress_blue = u[:, :rank] @ np.diag(S[:rank])@Vt[:rank, :]   # @ is used for matrix multiplication
    compress_blue = compress_blue.astype(int)
    #print(compress_blue)
    
    
    # for combining the three channels
    compressed_array = np.stack((compress_red,compress_green,compress_blue), axis=2)
    return compressed_array
    # for displaying the compressed image
    #plt.imshow(compressed_array)
    
    

# with my own multiplication function
def compress_color_b(rank,R,G,B):
    """_summary_ : This function is used to compress the image by using the rank of the image

    Args:
        rank (int): The rank of the image
        R (array ): The matrix of the image because of red channel  
        B (array): The matrix of the image because of blue channel
        G (array): The matrix of the image because of green channel
    """
    # calling the SVD function
    U, S, Vt = compress_matrix(rank,R)  # for red channel
    u,S,Vt = compress_matrix(rank,B)   # for blue channel
    u,s,Vt = compress_matrix(rank,G)   # for green channel
    
    a = matrix_product(U,S)
    compress_red = np.array(matrix_product(a,Vt)) 
    compress_red = compress_red.astype(int)
    #print(compress_red)
    b = matrix_product(u,S)
    compress_blue = np.array(matrix_product(b,Vt))  
    compress_blue = compress_blue.astype(int)
    #print(compress_blue)
    c = matrix_product(u,s)
    compress_green = np.array(matrix_product(c,Vt))  
    compress_green = compress_green.astype(int)
    #print(compress_green)
    
    # for combining the three channels
    compressed_array = np.stack((compress_red,compress_green,compress_blue), axis=2)
    return compressed_array
    # for displaying the compressed image
    # # plt.imshow(compressed_array)