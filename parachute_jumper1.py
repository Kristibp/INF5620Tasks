from vertical_motion import *
from numpy import *
from vertical_motion_force_plot import *

def input_parachute_jumper1():
    #Parachute jumper1 falling in air.

    d = 0.5 #m
    A = 0.9 #m^2
    V = 0.08 #m^3

    C_D = 1.2
    mu = 1.8*10**-5 #Pas
    rho = 0.79 #kg/m^3
    rho_b = 1003 #kg/m^3

    return C_D, d, mu, rho,rho_b, A, V

def run_parachute_jumper(T):
    v_0 = 0 #m/s
    g = 9.81 #m/s^2
    dt = 10**-1 #s
    force_quad = False
    v = vertical_motion(T, dt, v_0, g,input_terminal_v_test, force_quad)
    return v

#Was uncertain wether I should include output for this program, but as I found an unfixed bug in the program I included my testcase
v = run_parachute_jumper(10)
force_plot(v,9.81,10**-1,input_parachute_jumper1)
print 'v at T = 10:', v[-1]
