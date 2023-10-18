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