import matplotlib.pyplot as plt
from numpy import sqrt, log
from decimal import Decimal as Dc
from decimal import *
from math import ceil

getcontext().prec = 600

YSize = 9    # maximal allowable graph window width
XSize = 12   # maximal allowable graph window height
Draw_bubbles_in_scale = False  # if True, bubbles will be drawn in their real proportions
                               # if False, scale will be chosen by matplotlib
                               # try to experiment with it, that`s safe



MaxPower = 36     # determines what maximal power of theta will be taken into account
Margin = 10**-2   # maximal allowable error (Delta)
                  # if Delta > Margin, the calculation is discarded

Eta = Dc(1.2)     # Eta = Delta_z / R
L = Dc(0.7)       # L is Lambda
H = (Dc(1) - L) * Dc(0.999999999)



def zeta(z, fun, order):   # this one is very messy
    out = Dc(0)
    
    if(order == 2):
        for n in range(int((z + 1) / 2)):
            out += Dc(2) * fun(n) * fun(z - n)
            
        if(z % 2 == 0):  # if z is even
            out += fun(int(z / 2)) ** 2
        
    elif(order == 4):
        for a in range(int((z + 1) / 2)):
            f_z_a_2 = fun(int((z - a) / 2)) ** 2
            
            for b in range(int((a + 1) / 2)):
                f_a_b = fun(a - b)
                f_b = fun(b)
                
                for c in range(int((z - a + 1) / 2)):
                    out += Dc(8) * f_a_b * f_b * fun(c) * fun(z - a - c)
                
                if((z - a) % 2 == 0):
                    out += Dc(4) * f_a_b * f_b * f_z_a_2
            
            if(a % 2 == 0):
                f_a_2_2 = fun(int(a/2)) ** 2
                for c in range(int((z - a + 1) / 2)):
                    out += Dc(4) * f_a_2_2 * fun(c) * fun(z - a - c)
                    
                if((z - a) % 2 == 0):
                    out += Dc(2) * f_a_2_2 * f_z_a_2             
                    
        if(z % 2 == 0):  # if z is even
            if(z % 4 == 0):
                f_z4_2 = fun(int(z/4)) ** 2
                
            for b in range(int(z/2) + 1):
                f_b = fun(b)
                f_z2_b = fun(int(z/2) - b)
                
                for c in range(int((int(z/2) + 1) / 2)):
                    out += Dc(2) * f_z2_b * f_b * fun(c) * fun(int(z/2) - c)
                    
                if(z % 4 == 0):
                    out += f_z2_b * f_b * f_z4_2
        
    else:
        raise RuntimeError
            
    return out

def Str(n, l):
    if(n == 0):
        lg = 0
        f = 0.0
    elif(n == 1):
        lg = 0
        f = 1.0
    else:
        try:
            lg = log(abs(float(n))) / log(10)
            if(lg >= 0):
                lg = int(lg)
            else:
                lg = -ceil(-lg)
            f = float(n) / 10**lg
        except Exception:
            lg = 0
            f = 0
            print(n)
    
    out = ''
    out += (' ', '-')[f < 0]
    out += str(abs(round(f, l)))
    out += '0' * (2 + l - len(str(abs(round(f, l)))))
    out += 'e'
    out += ('+', '-')[lg < 0]
    out += str(abs(lg))
    out += ' ' * (3 - len(str(abs(lg))))
    
    return out

def Show_C(C):
    print('Taylor approximation for Y:')
    print('   Y ~ ')
    for z in range(len(C)):
        print('      ', end='')
        Print(C[z], 3)
        print('   theta ** ', end='')
        print(z)
        
def Print(n, l):
    print(Str(n, l), end='')

def calculate_C(l, h):
    C = [Dc(0)] * (MaxPower + 1)
    C[0] = Dc(1)
    C[1] = sqrt((Dc(1) / (l + h))**2 - Dc(1))
    A = Dc(-2) * C[1] * (l + h)**2
    B = l**2
    D = Dc(2) * l * h
    E = Dc(2) * l * h - Dc(1)
    F = (l + h)**2
    for z in range(1, MaxPower):
        k = tuple(C)
        
        S1 = Dc(0)
        for n in range(z):
            S1 += zeta(n, lambda i: Dc(i+1)*k[i+1], 2) * (B * zeta(z - n, lambda i: k[i], 4) + D * zeta(z - n, lambda i: k[i], 2))
        
        S2 = Dc(0)
        for n in range(int((z - 1) / 2)):
            S2 += Dc(2) * k[n + 2] * k[z - n] * Dc(n + 2) * Dc(z - n)
        if(z % 2 == 0):
            S2 += (k[int(z / 2) + 1] * Dc(int(z / 2) + 1)) ** 2
        S2 *= F
        
        c_new = (S1 + S2 + B * zeta(z, lambda i: k[i], 4) + E * zeta(z, lambda i: k[i], 2)) / A / Dc(z + 1)
        
        C[z + 1] = c_new
        
    return tuple(C)
    
def get_Y(theta, C):
    out = Dc(0)
    
    theta_pow_i = Dc(1)
    for i in range(len(C)):
        out += C[i] * theta_pow_i
        
        theta_pow_i *= theta
    
    return out

def get_Y_slope(theta, C):
    out = Dc(0)
    
    theta_pow_i = Dc(1)
    for i in range(len(C) - 1):
        out += (i + 1) * C[i + 1] * theta_pow_i
        
        theta_pow_i *= theta
    
    return out    

def draw_bubble(C, l, h, mu, eta, show):
    Xlist = []
    Zlist = []
    theta = Dc(0)
    dh = Dc(10**-3)
    y = get_Y(theta, C)
    y0 = get_Y(eta/2 * mu, C) / mu
    while(theta < eta/Dc(2) and y > 0):
        y = get_Y(theta * mu, C) / mu
        
        Xlist.append(float(y / y0))
        Zlist.append(float(theta))
        
        theta += dh
        
    errmax = get_errmax(C, l, h, eta * mu)
    
    if(Draw_bubbles_in_scale):
        scaleX = XSize  / (2 * max(Xlist))
        scaleY = YSize  / (max(Zlist) - min(Zlist))
        scale = min(scaleX, scaleY)
        X = scale * max(Xlist)
        Y = scale * (max(Zlist) - min(Zlist))
        
        fig = plt.gcf()
        fig.set_size_inches(X, Y)
        axs = fig.add_axes([ 0.1, 0.1, 0.8, 0.8])    
    else:
        axs = plt.gcf().add_axes([ 0.1, 0.1, 0.8, 0.8])
    
    axs.plot([0], [0])
    axs.plot(Xlist, Zlist, color='b')
    axs.plot(list(map(lambda Z: -Z, Xlist)), Zlist, color='b')
    axs.plot(Xlist, list(map(lambda Z: -Z, Zlist)), color='b')
    axs.plot(list(map(lambda Z: -Z, Xlist)), list(map(lambda Z: -Z, Zlist)), color='b')    
    plt.xlabel('r / R')
    plt.ylabel('Z / R')
    plt.title('Slice of a bubble\n' + r'$\Lambda = $' + str(round(l, 3)) + r'   $\eta = $' + str(round(eta, 3)) + r'   $\mu = $' + str(round(mu, 3)) + '   ' + r'$|\Delta|\eqslantless$ ' + Str(abs(errmax), 3))
    
    
    if(show):
        plt.show()   

def plot_Y(C, l, h, eta, show):
    Xlist = []
    Zlist = []
    theta = Dc(0)
    dh = Dc(10**-3)
    y = get_Y(theta, C)
    while(theta < eta/2):
        y = get_Y(theta, C)
        
        Xlist.append(float(theta))
        Zlist.append(float(y))
        
        theta += dh
        
    errmax = get_errmax(C, l, h, eta)
    
    plt.plot(Xlist, Zlist, color='b')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\Upsilon$')
    plt.title(r'Graph of $\Upsilon$' + '\n' + r'$\Lambda = $' + str(round(l, 3)) + r'   $\eta = $' + str(round(eta, 3)) + '   ' + r'$|\Delta|\eqslantless$ ' + Str(abs(errmax), 3))
    if(show):
        plt.show()   

def plot_error(C, l, h, eta):
    errlist = []
    tlist = []
    errmax = 0
    
    theta = 0
    dt = eta / 1000
    while(theta < eta/2):
        try:
            r = get_Y(theta, C)
            r_slope = get_Y_slope(theta, C)
            
            err = r / sqrt(1 + r_slope**2) - l * r**2 - h
            errlist.append(err)
            
            if(abs(err) > abs(errmax)):
                errmax = err
                
        except Exception:
            pass
        
        tlist.append(theta)
        theta += dt  
    
    plt.plot(tlist, errlist)
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\frac{\Upsilon}{\sqrt{1 + \dot{\Upsilon}^2}} - \Lambda\Upsilon^2 - H$')
    plt.title(r'Graph of deviation of $\Upsilon$' + '\n' + '$\Lambda$ = ' + str(round(l, 3)) + '   H = ' + str(round(h, 3)))    
    plt.show()
    
    print('Maximal error: ', end='')
    Print(errmax, 3)
    
def get_errmax(C, l, h, eta):
    errmax = 0
    
    theta = 0
    dt = eta / 200
    while(theta < eta/2):
        try:
            r = get_Y(theta, C)
            r_slope = get_Y_slope(theta, C)
            
            err = r / sqrt(1 + r_slope**2) - l * r**2 - h
            
            if(abs(err) > abs(errmax)):
                errmax = err
                
        except Exception as exc:
            print('Exception in get_errmax')
            input(exc)
        
        theta += dt  
    
    return errmax  
    
def boundary_dif(mu, theta, C):
    err = get_Y(theta * mu, C) - mu
    
    return err

def find_intervals(eta, C):
    points = []
    
    N = Dc(16)
    mu_min = Dc(0.1)
    mu_max = Dc(4)
    dmu = (mu_max - mu_min) / N
    mu = mu_min
    err = boundary_dif(mu, eta/Dc(2), C)
    last_sign = err / abs(err)
    while(mu < mu_max):
        mu += dmu
        
        err = boundary_dif(mu, eta/Dc(2), C)
        sign = err / abs(err)
        if(sign != last_sign):
            if(sign == 1):
                mu_n = mu - dmu
                mu_p = mu
            else:
                mu_n = mu
                mu_p = mu - dmu
                
            points.append((mu_n, mu_p))
            last_sign = sign
            print('    found an x-intercept between ' + Str(mu_n, 3) + ' and ' + Str(mu_p, 3))        
    print()
            
    return points
    
def calculate_mu(eta, C, l, h):
    out = []
    
    points = find_intervals(eta, C)
    
    for interval in points:
        mu_n = interval[0]
        mu_p = interval[1]
        
        print('    finding precise value of mu between ' + Str(mu_n, 3) + ' and ' + Str(mu_p, 3))
        for iteration in range(24):
            mid = (mu_n + mu_p) / Dc(2)
            
            err = float(boundary_dif(mid, eta/Dc(2), C))
            sign = err / abs(err)         
            
            if(sign == 1):
                mu_p = mid
            else:
                mu_n = mid
            
        if(abs(get_errmax(C, l, h, eta * mid)) < Margin):
            out.append(mid)
    
    return out

def draw_multiple_bubbles(C, l, h, mulist, eta):
    errmax = 0
    for mu in mulist:
        mu = Dc(mu)
        
        Xlist = []
        Zlist = []
        theta = Dc(0)
        dh = Dc(10**-3)
        y = get_Y(theta, C)
        y0 = get_Y(eta/2 * mu, C) / mu
        while(theta < eta/Dc(2) and y > 0):
            y = get_Y(theta * mu, C) / mu
            
            Xlist.append(float(y / y0))
            Zlist.append(float(theta))
            
            theta += dh
            
        errmax = max(errmax, abs(get_errmax(C, l, h, eta * mu)))
        
        plt.plot(Xlist, Zlist, color='b')
        plt.plot(list(map(lambda Z: -Z, Xlist)), Zlist, color='b')
        plt.plot(Xlist, list(map(lambda Z: -Z, Zlist)), color='b')
        plt.plot(list(map(lambda Z: -Z, Xlist)), list(map(lambda Z: -Z, Zlist)), color='b') 
        
    possible_mu = ''
    for mu in mulist[:-1]:
        possible_mu += str(round(mu, 3))
        possible_mu += ' or '
    possible_mu += str(round(mulist[-1], 3))
    
    plt.plot([0], [0])
    plt.xlabel('r / R')
    plt.ylabel('Z / R')
    plt.title('Slice of a bubble\n' + r'$\Lambda = $' + str(round(l, 3)) + r'   $\eta = $' + str(round(eta, 3)) + r'   $\mu = $' + possible_mu + '   ' + r'$|\Delta|\eqslantless$ ' + Str(abs(errmax), 3))
    
    plt.show()           
    
    
    
def main():
    print('Calculating C')
    C = calculate_C(L, H)
    
    Show_C(C)
    
    mulist = calculate_mu(Eta, C, L, H)
    if(len(mulist) == 0):
        print('Found no solutions!')
    elif(len(mulist) == 1):
        print('Found one solution')
    else:
        print('Found ' + str(len(mulist)) + ' solutions')
        
    for mu in mulist:    
        draw_bubble(C, L, H, mu, Eta, True)
        plot_Y(C, L, H, Eta, True)
        plot_error(C, L, H, Eta)
    draw_multiple_bubbles(C, L, H, mulist, Eta)
    
    print()
    print('Program finished')
    input()

try:
    main()
except Exception as exc:
    input(exc)
