from numpy import *
import nose.tools as nt

#Note: I think there is an error implementation in this program. The jumper in task 26 is falling way to slow... and in the wrong direction too. I couldn't get the Stokes method to pass all the nose tests either. I am however running out of time, so I hand in the program anyway. Thought it would be better to show I've done something instead of handing in nothing.

def compute_stokes_a(d,mu,rho_b,V):
    a = 3*pi*d*mu/float(rho_b*V)
    return a

def compute_b(g,rho, rho_b):
    b = g*(rho/float(rho_b) - 1)
    return b

def compute_quad_a(rho, C_D, A, rho_b, V):
    a = 0.5*C_D*rho*A/float(rho_b*V)
    return a

def compute_Re(rho, d,v,mu):
    return rho*d*abs(v)/float(mu)

def step(C_D, d, mu, rho,rho_b, A, V, g, v_n,dt,force_quad):

    def stokes_step():
        a  = compute_stokes_a(d,mu,rho_b,V)
        v = ((1 - dt*a*0.5)*v_n + dt*b)/float(1 + dt*a*0.5)
        return v

    def quad_step():
        a = compute_quad_a(rho, C_D, A, rho_b, V)
        v = (v_n + dt*b)/float(1 + dt*a*abs(v_n))
        return v
    
    Re = compute_Re(rho, d,v_n,mu)
    b = compute_b(g,rho, rho_b)
    if Re > 1 or force_quad: v = quad_step()
    else: v = stokes_step()
    return v

def read_command_line():
    import sys
    if len(sys.argv) < 2:
        sys.exit(1)
    args = []
    for i in range(1,len(sys.argv)):
        args.append(float(sys.argv[i]))
    return args

def input_heavy_sphere():
    #Heavy sphere of radius r = 5 falling in air.
    
    r = 5 #m
    d = 2*r #m
    A = pi*r**2 #m^2
    V = 4*pi*r**3/3 #m^3

    C_D = 0.45
    mu = 1.8*10**-5 #Pas
    rho = 0 #kg/m^3
    rho_b = 100000 #kg/m^3

    return C_D, d, mu, rho,rho_b, A, V

def input_terminal_v_test():
    #Heavy sphere of radius r = 5 falling in air.
    
    r = 5 #m
    d = 2*r #m
    A = pi*r**2 #m^2
    V = 4*pi*r**3/3 #m^3

    C_D = 0.45
    mu = 1.8*10**-5 #Pas
    rho = 10 #kg/m^3
    rho_b = 9 #kg/m^3

    return C_D, d, mu, rho,rho_b, A, V

def test_linear(): 
    #Note: In this test I don't test Stokes implementation. I don't see why F_d^S would be 0 if rho = 0
    v_0 = 0 #m/s
    T = 10 #s
    g = 9.81 #m/s^2
    dt = 1 #s
    force_quad = True
    v = vertical_motion(T, dt, v_0, g,input_heavy_sphere,force_quad)
    v_T = v_0 - g*T
    nt.assert_almost_equal(v[-1] - v_T, 0, 14)

def test_terminal_v():
    v_0 = 10 #m/s
    T = 100 #s, a value I found through testing was high enough for the solution to clearly converge.
    g = 9.81 #m/s^2
    dt = 10**-1 #s
    force_quad = False
    v = vertical_motion(T, dt, v_0, g,input_terminal_v_test, force_quad)

    C_D, d, mu, rho,rho_b, A, V = input_terminal_v_test()
    v_T_1 = compute_b(g,rho, rho_b)*compute_stokes_a(d,mu,rho_b,V)
    v_T_2 = sqrt(compute_b(g,rho, rho_b)/compute_quad_a(rho, C_D, A, rho_b, V))
    diff = min(abs(v[-1]-v_T_1), abs(v[-1])-abs(v_T_2)) #Finding the best approximation of the analytic values
    nt.assert_almost_equal(diff, 0, 12)

def vertical_motion(T, dt, v_0, g, f, force_quad):
    C_D, d, mu, rho,rho_b, A, V = f()
    N = int(T/dt)
    v = zeros(N+1)
    v[0] = v_0
    for i in range(N):
        v[i+1] = step(C_D, d, mu, rho,rho_b, A, V, g, v[i], dt, force_quad)
    return v
