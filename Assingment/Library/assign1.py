#### Assingment 1 

import matplotlib.pyplot as plt 


## 1. Sum of first N odd numbers
def odd (N):
    # Checking the input
    if N<0:
        raise ValueError("N should be a natural number")
    
    series = []  #create an empty list
    for i in range(N):
        series.append(2*i+1)
    return series

### Factorial of N
def fac(n):
    x = 1
    for i in range(1,n+1):
        x = x*i
    return x




## 2. Calculate the sum of N terms of an AP series
def AP(N, d: float = 1.5, a: float = 1):
    """Finds the sum of an Arithmatic Progression.

    Args:
        N (int, optional): N is the number of terms to return. 
        d (float, optional): common difference of the AP series. Defaults to 2.
        a (float, optional): first term of the AP seies. Defaults to 1.

    Raises:
        ValueError: when N is less than 0

    Returns:
        list: AP seies of N terms 
    """                                 
    
    # Checking the input
    if N<0:
        raise ValueError("N should be a natural number")
    
    series = []  #create an empty list
    for i in range(N):
        series.append(a+i*d)
    return series


# Calculate the sum of N terms of a GP series
def GP(N, r: float = 0.5, a: float = 1):
    """Finds the sum of an Arithmatic Progression.

    Args:
        N (int, optional): N is the number of terms to return. 
        r (float, optional): common ratio of the GP series. Defaults to 0.5.
        a (float, optional): first term of the GP seies. Defaults to 1.

    Raises:
        ValueError: when N is less than 0

    Returns:
        list: GP seies of N terms 
    """                                 
    
    # Checking the input
    if N<0:
        raise ValueError("N should be a natural number")
    
    series = []  #create an empty list
    for i in range(N):
        series.append(a*r**i)
    return series


# Calculate the sum of N terms of a HP series
def HP(N, d: float = 1.5, a: float = 1):
    """Finds the sum of an Arithmatic Progression.

    Args:
        N (int, optional): N is the number of terms to return. 
        d (float, optional): common difference of the HP series. Defaults to 0.5.
        a (float, optional): first term of the HP seies. Defaults to 1.

    Raises:
        ValueError: when N is less than 0

    Returns:
        list: GP seies of N terms 
    """                                 
    
    # Checking the input
    if N<0:
        raise ValueError("N should be a natural number")
    
    series = []  #create an empty list
    for i in range(N):
        series.append(1/(a+i*d))
    return series



## 3. Sum of power series
def ser(N):
 if N >0 :
    # create an empty list
    x_list = []  
    y_list = []
    s = 0
    for i in range(1,N+1):
        x_list.append(i)
        s = s + ((-1)**(i+1)/2**i)
        y_list.append(s)
 else :
   raise ValueError ("N should be a natural number") 
 print(f"sum of all the terms is ",round(sum(y_list[-2:]),4))

 plt.plot (x_list,y_list)
 plt.xlabel ("N")
 plt.ylabel ("Sum")
 plt.title ("Sum vs N graph")
 


## 4. Matrix Multiplication
# A.B
def mat1(A,B): 
 r = [[0,0,0],
    [0,0,0],
    [0,0,0],]
 for i in range(len(A)):  # iterating by row of A
    # iterating by column by B
    for j in range(len(B[0])):   
    # iterating by row of B
        for k in range(len(B)):  
            r[i][j] += A[i][k]*B[k][j]
            
 for m in r:
   print(m)

# C.D
def mat2(C,D):  
 r = [[0]]

 for i in range(len(C)):        # iterating by row of C
    # element wise multiplication 
    r[0][0] += (C[i][0])*(D[i][0])    
    
 print(r)  
 
# B.C
def mat3(B,C): 
 r = [[0],
    [0],
    [0],]
 for i in range(len(B)):               # iterating by row of B
    for j in range(len(C[0])):             # iterating by column by C
        for k in range(len(C)):                    # iterating by row of C
            r[i][j]+= B[i][k]*C[k][j]
 for m in r:
    print(m)
 
 
 
 
 ## Class for complex numbers
class myComplex:
    def __init__(self,a,b,x,y) :       # a & x are real part and b & y are imaginary part
      # instance variables
      self.a = a
      self.b = b
      self.x = x
      self.y = y
      self.mod_1 = (self.a**2 + self.b**2)**(1/2)
      self.mod_2 = (self.x**2 + self.y**2)**(1/2)
      
    # Add complex numbers
    def Add(self):
        self.sum_real = self.a + self.x
        self.sum_imaginary = self.b + self.y
        return (print(f"sum of complex number is {self.sum_real} + {self.sum_imaginary}i"))
    
    # product of complex numbers 
    def Multiply(self):
        self.product_real = self.a*self.x - self.b + self.y
        self.product_imaginary = self.b + self.x + self.a*self.y
        return (print(f"product of complex number is {self.product_real} + {self.product_imaginary}i"))  
    
    # modulus of complex numbers
    def Modulus(self):
        return (print(f"modulus of complex number {self.a} + {self.b}i is {self.mod_1} \nmodulus of complex number {self.x} + {self.y}i is {self.mod_2}"))
