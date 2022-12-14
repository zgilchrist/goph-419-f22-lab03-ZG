import numpy as np

def ode_freefall_euler(g0, dg_dz, cd_star, H, dt):
    """Solves a second-order ode using Euler's Method(RK1)
    
    Parameters
    ----------
    g0: float
        Gravity constant g_0

    dg_dz: float
        Free-air gradient of gravitational acceleration g'

    cd_star: float
        Mass normalized drag coefficient c_D* 
    
    H: float
        Total drop height H
    
    dt: float
        Time step delta t 
    Returns
    -------
    t: np.ndarray
        Array of shape(n) that contains time t values
    z: np.ndarray
        Array of shape(n) that contains distance z(t) values
    v: np.ndarray
        Array of shape(n) that contains velocity v(t) values 
    """

    #find N using N=(b-a)/h to find difference 
    a = 0
    b = H
    N = int(H/dt)
    
    t = np.zeros(N)
    z = np.zeros(N)
    v = np.zeros(N)
    a = np.zeros(N)

    t[0] = 0
    z[0] = 0
    v[0] = 0
    a[0] = 0

    #z[-1] = H
    
    for i in np.arange(1,N):
        t[i] = dt * i
        a[i] = a[i-1] + dt * (g0 + dg_dz * z[i-1])
        v[i] = v[i-1] + dt * (g0 + dg_dz * z[i-1] - a[i-1]) / cd_star
        z[i] = z[i-1] + dt * -1 * (g0 - a[i-1] - cd_star * v[i-1]) / cd_star
        
    return t, z, v

def ode_freefall_rk4(g0, dg_dz, cd_star, H, dt):
    """Solves a second-order ode using the classical RK4 Method from the set of Runge-Kutta methods 
    
    Parameters
    ----------
    g0: float
        Gravity constant g_0

    dg_dz: float
        Free-air gradient of gravitational acceleration g'

    cd_star: float
        Mass normalized drag coefficient c_D*
    
    H: float
        Total drop height H
    
    dt: float
        Time step delta t 
    Returns
    -------
    t: np.ndarray
        Array of shape(n) that contains time t values
    z: np.ndarray
        Array of shape(n) that contains distance z(t) values
    v: np.ndarray
        Array of shape(n) that contains velocity v(t) values 
    """
    t = np.zeros()
    z = np.zeros()
    v = np.zeros()

    t[0] = 0
    z[0] = 0
    v[0] = 0
    return t, z, v