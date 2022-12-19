"""Driver script for Assignment #3"""
#F22 GOPH 419 Brandon Karchewski
#Zachary Gilchrist

from time import perf_counter_ns
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import ode_solver

#constants for the equations
#globally accessible
g0 = 9.811636
dg_dz = 3.086e-6
mu_0 = 1.827e-5
diameter = 0.015
radius = diameter/2
rho_steel = 7800
vol_bb = 4/3*np.pi*radius**3
mass_bb = vol_bb*rho_steel
cd = 6*np.pi*mu_0*radius
cd_star = cd/mass_bb

def main():
    """Main function"""

    #sample value drop heights
    H = np.array([10,20,40])
    
    #time step in seconds
    dt = np.array([1e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.25,0.5,1,2,3,10])
    
    ##############################################

    #q1 t* dt plots

    plt.figure()
    plt.subplot(311)
    plt.loglog(dt,get_t_star(dt,H[0],'rk1'),'r-o',label = 'RK1')
    plt.loglog(dt,get_t_star(dt,H[0],'rk4'),'b-o',label = 'RK4')
    plt.title('H = 10')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")
    
    plt.subplot(312)
    plt.loglog(dt,get_t_star(dt,H[1],'rk1'),'r-o',label = 'RK1')
    plt.loglog(dt,get_t_star(dt,H[1],'rk4'),'b-o',label = 'RK4')
    plt.title('H = 20')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")
    
    plt.subplot(313)
    plt.loglog(dt,get_t_star(dt,H[2],'rk1'),'r-o',label = 'RK1')
    plt.loglog(dt,get_t_star(dt,H[2],'rk4'),'b-o',label = 'RK4')
    plt.title('H = 40')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")

    ##############################################

    #q2 performance time dt plots
    
    plt.figure()
    plt.loglog(dt,get_t_star(dt,H[0],'rk4',True),'r-o',label ='H=10')
    plt.loglog(dt,get_t_star(dt,H[1],'rk4',True),'b-o',label ='H=20')
    plt.loglog(dt,get_t_star(dt,H[2],'rk4',True),'g-o',label ='H=40')
    plt.title('RK4 for three drop heights')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Simulation Time (ns)")

    ##############################################

    #q3 sensitivity analysis plots

    
    plt.show()


def get_t_star(dt,H,method,get_time = False):
    if get_time == True:
        t_star = np.zeros(0)
        time = np.zeros(0)
        if method == 'rk1':
            for i in dt:
                start = perf_counter_ns()
                t,z,v = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H, i)
                t_star = np.append(t_star,t[-1])
                end = perf_counter_ns()
                time = np.append(time,end-start)
        if method == 'rk4':
            for i in dt:
                start = perf_counter_ns()
                t,z,v = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H, i)
                t_star = np.append(t_star,t[-1])
                end = perf_counter_ns()
                time = np.append(time,end-start)
        return time
    
    elif get_time == False:
        t_star = np.zeros(0)
        if method == 'rk1':
            for i in dt:
                t,z,v = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H, i)
                t_star = np.append(t_star,t[-1])
        if method == 'rk4':
            for i in dt:
                t,z,v = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H, i)
                t_star = np.append(t_star,t[-1])
        return t_star


if __name__ == "__main__":
        main()