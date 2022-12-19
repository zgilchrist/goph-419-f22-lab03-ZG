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
    
    #N = int(H/dt)
    N = 1
    t = np.zeros(N)
    z = np.zeros(N)
    v = np.zeros(N)
    
    t[0] = 0
    z[0] = 0
    v[0] = 0
    
    i=1
    while z[-1] < H:
        t = np.append(t,t[i-1] + dt)
        #print('Value of t[i] #1: ',t[i])
        v = np.append(v,v[i-1] + dt * (g0 + dg_dz * z[i-1] - cd_star * v[i-1]))
        #print('Value of v[i] #1: ', v[i])
        z = np.append(z,z[i-1] + dt * v[i-1])
        #print('Value of z[i] #1: ',z[i])
        

        if z[i] >= H:
            z[i] = H
            #print('Value of z[i] #2: ',z[i])
            dt = (z[i] - z[i-1])/v[i-1]
            t = np.append(t,t[i-1] + dt)
            #print('Value of t[i] #2: ', t[i])
            v = np.append(v,v[i-1] + dt * (g0 + dg_dz * z[i-1] - cd_star * v[i-1]))
            #break
            
        #print('Value of i: ',i)
        i+=1
        
        
    #print("t is: ",t)
    #print('Value of t[-1]: ',t[-1])
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
    #find N using N=(b-a)/h to find difference 
    
    #N = int(H/dt)
    N = 1
    t = np.zeros(N)
    z = np.zeros(N)
    v = np.zeros(N)
    
    t[0] = 0
    z[0] = 0
    v[0] = 0
    
    i=1
    while z[-1] < H:
        t = np.append(t,t[i-1] + dt)
        #print('Value of t[i] #1: ',t[i])
        k1_v = (g0 + dg_dz * z[i-1] - cd_star * v[i-1])
        k1_z = v[i-1]
        k2_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + 0.5 * k1_v * dt))
        k2_z = (0.5 * k1_z * dt) + v[i-1]
        k3_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + 0.5 * k2_v * dt))
        k3_z = (0.5 * k2_z * dt) + v[i-1]
        k4_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + k3_v * dt))
        k4_z = (k3_z * dt) + v[i-1]
        
        rk4_v = v[i-1] + (k1_v+2*k2_v+2*k3_v+k4_v)*dt/6
        rk4_z = z[i-1] + (k1_z+2*k2_z+2*k3_z+k4_z)*dt/6

        v = np.append(v,rk4_v)
        #print('Value of v[i] #1: ', v[i])
        z = np.append(z,rk4_z)
        #print('Value of z[i] #1: ',z[i])
        

        if z[i] >= H:
            z[i] = H
            #print('Value of z[i] #2: ',z[i])
            dt = (z[i] - z[i-1])/v[i-1]
            t = np.append(t,t[i-1] + dt)
            #print('Value of t[i] #2: ', t[i])

            k1_v = (g0 + dg_dz * z[i-1] - cd_star * v[i-1])
            k2_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + 0.5 * k1_v * dt))
            k3_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + 0.5 * k2_v * dt))
            k4_v = (g0 + dg_dz * z[i-1] - cd_star * (v[i-1] + k3_v * dt))
            
            rk4_v = v[i-1] + (k1_v+2*k2_v+2*k3_v+k4_v)*dt/6
        
            v = np.append(v,rk4_v)
            
        #print('Value of i: ',i)
        i+=1
        
        
    #print("t is: ",t)
    #print('Value of t[-1]: ',t[-1])
    return t, z, v
