from vertical_motion import *
from numpy import *
import matplotlib.pyplot as plt

def compute_forces(v,g,f):
    C_D, d, mu, rho,rho_b, A, V = f()

    F_g = -rho_b*V*g
    F_b = rho*V*g
    F_Ds = 3*pi*d*mu*v
    F_Dq = -0.5*C_D*rho*A*abs(v)*v
    return F_g, F_b, F_Ds, F_Dq

def test_plot_forces():
    v_0 = 10 #m/s
    T = 10 #s
    g = 9.81 #m/s^2
    dt = 10**-3 #s
    force_quad = False
    v = vertical_motion(T, dt, v_0, g,input_terminal_v_test, force_quad)
    force_plot(v,g,dt,input_terminal_v_test)

def force_plot(v,g,dt,input):
    F_g, F_b, F_D,t = zeros(len(v)), zeros(len(v)), zeros(len(v)),zeros(len(v))
    C_D, d, mu, rho,rho_b, A, V = input()

    for i in range(len(v)):
       F_g[i], F_b[i], F_Ds, F_Dq = compute_forces(v[i],g,input)
       if compute_Re(rho, d,v[i],mu) > 1: F_D[i] = F_Dq
       else: F_D[i] = F_Ds
       t[i] = dt*i

    plt.hold('off')
    pg = plt.plot(t, F_g, 'r-')
    plt.hold('on')
    pb = plt.plot(t, F_b, 'g-') 
    pD = plt.plot(t, F_D, 'b-')
    plt.legend(['F_g', 'F_b', 'F_D'])
    plt.show()
    plt.hold('off')

if __name__ == '__main__':
    test_plot_forces()
