
### Solving the differential equation using the Runge-Kutta method

def rk4(d2ydx2, dydx, x0, y0, z0, xf, h):
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