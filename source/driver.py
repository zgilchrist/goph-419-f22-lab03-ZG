"""Driver script for Assignment #3"""
#F22 GOPH 419 Brandon Karchewski
#Zachary Gilchrist

import numpy as np
import matplotlib.pyplot as plt

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
    dt = 0.001

    
    #for i in H:
    t1,z1,v1 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[0], dt)
    t2,z2,v2 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[1], dt)
    t3,z3,v3 = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H[2], dt)
    
    
    #ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, i, dt)

    plt.figure()
    plt.subplot(311)
    plt.plot(t1,z1,'r',t1,v1,'b')
    plt.subplot(312)
    plt.plot(t2,z2,'r',t2,v2,'b')
    plt.subplot(313)
    plt.plot(t3,z3,'r',t3,v3,'b')
    #plt.plot(H,ys_rk4)
    plt.show()







if __name__ == "__main__":
        main()