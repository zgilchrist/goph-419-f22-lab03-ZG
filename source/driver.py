"""Driver script for Assignment #3"""
#F22 GOPH 419 Brandon Karchewski
#Zachary Gilchrist

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import ode_solver

def main():
    """Main function"""

    #constants for the equations
    g0 = 9.811636
    dg_dz = 0.3086e-6
    mu_0 = 1.827e-5
    diameter = 0.015
    radius = diameter/2
    rho_steel = 7800
    vol_bb = 4/3*np.pi*radius**3
    mass_bb = vol_bb*rho_steel
    cd = 6*np.pi*mu_0*radius
    cd_star = cd/mass_bb
    H = np.array([10,20,40])
    
    #time step in seconds
    dt = 0.1

    
    #for i in H:
    t1,z1,v1 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[0], dt)
    t2,z2,v2 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[1], dt)
    t3,z3,v3 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[2], dt)
    
    
   #t4,z4,v4 = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H[0], dt)
    #t5,z5,v5 = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H[1], dt)
    #t6,z6,v6 = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H[2], dt)

    #rk1 plots
    plt.figure()
    plt.subplot(321)
    plt.plot(t1,z1,'r',t1,v1,'b')
    plt.title('RK1, H = 10')
    plt.subplot(323)
    plt.plot(t2,z2,'r',t2,v2,'b',label='H=20')
    plt.title('RK1, H = 20')
    plt.subplot(325)
    plt.plot(t3,z3,'r',t3,v3,'b',label='H=20')
    plt.title('RK1, H = 40')

    #rk4 plots
    plt.subplot(322)
    #plt.plot(t4,z4,'r',t4,v4,'b')
    plt.title('RK4, H = 10')
    plt.subplot(324)
    #plt.plot(t5,z5,'r',t5,v5,'b')
    plt.title('RK4, H = 20')
    plt.subplot(326)
    #plt.plot(t6,z6,'r',t6,v6,'b')
    plt.title('RK4, H = 40')
    
    plt.show()







if __name__ == "__main__":
        main()