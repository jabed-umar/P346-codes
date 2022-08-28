  
import matplotlib.pyplot as plt

# Random Number Generator by LCG method
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
    rand=[]                # This is the list of random numbers
    # 100 random numbers
    if k==0:
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


# volume of a sphere by random number
def vol(n):
    """This function calculates the volume of a sphere using random number
    Args:
        n (integer): This is the number of random numbers in the list
    """
    x = My_random(10,n,0)
    y = My_random(5,n,0)
    z = My_random(3,n,0)
    inside = 0
    for i in range(n):
        if x[i]**2 + y[i]**2 + z[i]**2 <= 1:
            inside += 1
    vol_s = inside/n
    print(vol_s)
    return vol_s
    # print(len(x))
    
# 2d random walk 
def steps(n):
    """This function calculates the number of steps in a 2d random walk
    Args:
        n (integer): This is the number of random numbers in the list
    """
    rms_0 = 0
    x = My_random(10,n,1)
    y = My_random(5,n,1)
    plt.plot(x,y)
    plt.show()
    for i in range(1,n):
        rms_0 += ((x[i] - x[i-1])**2 + (y[i] - y[i-1])**2)
    print("RMS of step" ,n, "is: ",(rms_0/n)**0.5)
    dis = ((x[n-1]**2 + y[n-1]**2)**0.5)
    print("Displacement of step" ,n, "is: ",dis)
    com = rms_0/dis
    print("RMS of step" ,n, "is",com,"times the Displacement")