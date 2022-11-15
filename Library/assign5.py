
### integration by mid point method 
def midpoint(f, a, b, n):
    h = (b-a)/n
    s = 0
    for i in range(n):
        s += f(a + h*(i+0.5))
    # print ("For n =", n)
    return s*h


### integration by trapezoidal method
def trapezoidal(f, a, b, n):
    h = (b-a)/n
    s = 0.5*(f(a) + f(b))
    for i in range(1, n):
        s += f(a + h*i)
    print ("For n =", n)
    return s*h

#3)Simpson Methode
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

# integration by Monte Carlo method
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
