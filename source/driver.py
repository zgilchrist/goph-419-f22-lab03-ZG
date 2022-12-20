"""Driver script for Assignment #3"""
#F22 GOPH 419 Brandon Karchewski
#Zachary Gilchrist

from time import perf_counter_ns
import numpy as np
import matplotlib.pyplot as plt

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
    plt.title('Total drop time vs. time step at H = 10 in loglog space')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")
    
    plt.subplot(312)
    plt.loglog(dt,get_t_star(dt,H[1],'rk1'),'r-o',label = 'RK1')
    plt.loglog(dt,get_t_star(dt,H[1],'rk4'),'b-o',label = 'RK4')
    plt.title('Total drop time vs. time step at H = 20 in loglog space')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")
    
    plt.subplot(313)
    plt.loglog(dt,get_t_star(dt,H[2],'rk1'),'r-o',label = 'RK1')
    plt.loglog(dt,get_t_star(dt,H[2],'rk4'),'b-o',label = 'RK4')
    plt.title('Total drop time vs. time step at H = 40 in loglog space')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Total Drop Time (s)")
    plt.savefig('t_star_dt.png')

    plt.figure()
    plt.loglog(dt[1:],get_approx_re(get_t_star(dt,H[1],'rk1')),'r-o',label = 'RK1')
    plt.loglog(dt[1:],get_approx_re(get_t_star(dt,H[1],'rk4')),'b-o',label = 'RK4')
    plt.title('Appoximate Relative Error vs Time Step for i-1 elements at H=20')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Approximate Relative Error (%)")
    plt.savefig('approx_re.png')

    ##############################################

    #q2 performance time dt plots
    
    plt.figure()
    plt.loglog(dt,get_t_star(dt,H[0],'rk4',get_time = True),'r-o',label ='H=10')
    plt.loglog(dt,get_t_star(dt,H[1],'rk4',get_time = True),'b-o',label ='H=20')
    plt.loglog(dt,get_t_star(dt,H[2],'rk4',get_time = True),'g-o',label ='H=40')
    plt.title('RK4 for three drop heights in loglog space')
    plt.legend()
    plt.xlabel("Time Step (s)")
    plt.ylabel("Simulation Time (ns)")
    plt.savefig('simulation_time.png')
    ##############################################

    #q3 sensitivity analysis plots

    #g0
    plt.figure()
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.01),'r-o',label ='1%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.05),'b-o',label ='5%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.1),'g-o',label ='10%')
    plt.title('Percent Pertubation in g0 with RK4 at H=20')
    plt.legend()
    plt.xlabel("Time Step (s) in log space")
    plt.ylabel("Sensitivity of drop time (ms)")
    plt.xscale('log')
    plt.savefig('g0_sensitivity.png')

    #dg_dz
    plt.figure()
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',dg_dz_perturb=0.01),'r-o',label ='1%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',dg_dz_perturb=0.05),'b-o',label ='5%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',dg_dz_perturb=0.1),'g-o',label ='10%')
    plt.title("Percent Pertubation in g' with RK4 at H=20")
    plt.legend()
    plt.xlabel("Time Step (s) in log space")
    plt.ylabel("Sensitivity of drop time (ms)")
    plt.xscale('log')
    plt.savefig('dg_dz_sensitivity.png')

    #cd_star
    plt.figure()
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',cd_star_perturb=0.01),'r-o',label ='1%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',cd_star_perturb=0.05),'b-o',label ='5%')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',cd_star_perturb=0.1),'g-o',label ='10%')
    plt.title('Percent Pertubation in cd* with RK4 at H=20')
    plt.legend()
    plt.xlabel("Time Step (s) in log space")
    plt.ylabel("Sensitivity of drop time (ms)")
    plt.xscale('log')
    plt.savefig('cd_star_sensitivity.png')

    #effect of drop height on sensitivity
    plt.figure()
    plt.plot(dt,get_sensitivity(dt,H[0],'rk4',g0_perturb=0.01),'r-o',label ='1%, H=10')
    plt.plot(dt,get_sensitivity(dt,H[0],'rk4',g0_perturb=0.05),'b-o',label ='5%, H=10')
    plt.plot(dt,get_sensitivity(dt,H[0],'rk4',g0_perturb=0.1),'g-o',label ='10%, H=10')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.01),'c-o',label ='1%, H=20')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.05),'m-o',label ='5%, H=20')
    plt.plot(dt,get_sensitivity(dt,H[1],'rk4',g0_perturb=0.1),'y-o',label ='10%, H=20')
    plt.plot(dt,get_sensitivity(dt,H[2],'rk4',g0_perturb=0.01),'k-o',label ='1%, H=40')
    plt.plot(dt,get_sensitivity(dt,H[2],'rk4',g0_perturb=0.05),'b-d',label ='5%, H=40')
    plt.plot(dt,get_sensitivity(dt,H[2],'rk4',g0_perturb=0.1),'g-d',label ='10%, H=40')
    plt.title('Percent Pertubation in g0 with RK4 different drop heights')
    plt.legend()
    plt.xlabel("Time Step (s) in log space")
    plt.ylabel("Sensitivity of drop time (ms)")
    plt.xscale('log')
    plt.savefig('drop_height_sensitivity.png')
    
    plt.show()


def get_t_star(dt,H,method,get_time = False,):
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

#calculates error based on a given re
def get_approx_re(t_star):
    approx_re = np.zeros(0)
    i=1
    while i < len(t_star):
        approx_re = np.append(approx_re,100*(abs(t_star[i]-t_star[i-1])/t_star[i-1]))
        i+=1
    return approx_re

#method based off t_star that allows the user to modify the static parameters
def get_sensitivity(dt,H,method,g0_perturb = 0,dg_dz_perturb = 0,cd_star_perturb = 0):
    t_star = np.zeros(0)
    t_star_perturb = np.zeros(0)
    if method == 'rk1':
        for i in dt:
            t,z,v = ode_solver.ode_freefall_euler(g0, dg_dz, cd_star, H, i)
            t_star = np.append(t_star,t[-1])
            t,z,v = ode_solver.ode_freefall_euler(g0 + (g0*g0_perturb), dg_dz + (dg_dz*dg_dz_perturb), cd_star + (cd_star*cd_star_perturb), H, i)
            t_star_perturb = np.append(t_star_perturb,t[-1])
    if method == 'rk4':
        for i in dt:
            t,z,v = ode_solver.ode_freefall_rk4(g0, dg_dz, cd_star, H, i)
            t_star = np.append(t_star,t[-1])
            t,z,v = ode_solver.ode_freefall_euler(g0 + (g0*g0_perturb), dg_dz + (dg_dz*dg_dz_perturb), cd_star + (cd_star*cd_star_perturb), H, i)
            t_star_perturb = np.append(t_star_perturb,t[-1])
    #gets difference between original time and pertubated time in ms
    return np.subtract(t_star_perturb,t_star)*1000

if __name__ == "__main__":
        main()