"""
TODO
    Add sample code to run functions (savefig)
    
    waveguide class...
    file for solvers
    file for fields/confinement
    file for angles, film thicknesses
    file for loss
    utils file
"""

from numpy import *
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import *
# %matplotlib inline

"""
------------------------------------------------------------------------------------------------------------
INITIAL PARAMETERS
------------------------------------------------------------------------------------------------------------
"""
si = 3.5
sio2 = 1.45
air = 1
duty_cycle = 0.4

def refrac1(duty_cycle):
    """Estimates the effective refractive index of the film region of a swg waveguide

    PARAMETERS
    ----------
    duty_cycle : number
        subwavelength grating duty cycle
        0 ≤ duty_cycle ≤ 1

    RETURNS
    -------
    index : number
        refractive index

    """
    return si * duty_cycle + (1 - duty_cycle) * sio2

def refrac2(duty_cycle):
    """Estimates the effective refractive index of the film region of a swg waveguide

    PARAMETERS
    ----------
    duty_cycle : number
        subwavelength grating duty cycle
        0 ≤ duty_cycle ≤ 1
        
    RETURNS
    -------
    index : number
        refractive index

    """
    return sqrt(si**2 * duty_cycle + (1 - duty_cycle) * sio2**2)

def refrac3(duty_cycle):
    """Estimates the effective refractive index of the film region of a swg waveguide

    PARAMETERS
    ----------
    duty_cycle : number
        subwavelength grating duty cycle
        0 ≤ duty_cycle ≤ 1
        
    RETURNS
    -------
    index : number
        refractive index

    """
    return (si**-2 * duty_cycle + (1 - duty_cycle) * sio2**-2)**-0.5

def recover_duty_cycle(n):
    """Estimates the duty cycle of an SWG with the given effective index

    PARAMETERS
    ----------
    n : number
        effective refractive index

    RETURNS
    -------
    duty cycle : number


    NOTES
    -----
    • If it is not possible to create a polySi-SiO2 waveguide with this refractive index, it returns None

    """
    if n < sio2 or n > si:
        return

    return (n**2 - ns**2) / (nf**2 - ns**2)


# ASSIGNING REFRACTIVE INDEX OF THE SWG
nf = refrac2(duty_cycle)
nc = sio2
ns = sio2

n1 = nf
n2 = nc
n4 = ns
n3 = nc
n5 = ns

h = 200e-9 # height, STANDARD (FIXED) FOR BULK CMOS PROCESS
w = 450e-9 # width

default_h = h
default_w = w
fixed_h = h
fixed_w = w

"""
 --------------
|    | n2 |    |    n2: x (nc)
 --------------
| n5 | n1 | n3 |    n5: y (ns), n1: nf, n3: y (nc)
 --------------
|    | n4 |    |    n4: x (ns)
 --------------

 y(h)    z
    |   /
    |  /
    | /
    • ----- x(w)

"""

wavelength = 1550e-9
c_light = 3e8
f = c_light / wavelength
omega = f * 2 * pi # f = 2πw
k_lamb = lambda lamb: 2 * pi / lamb # evaluate k (wavenumber) for a given lambda
k = 2 * pi / wavelength
# beta = k * N
a = (ns**2 - nc**2) / (nf**2 - ns**2)
V = lambda lamb: k_lamb(lamb) * w * sqrt(nf**2 - ns**2) # evaluate V for a given lambda
V_nf = lambda lamb, nf: k_lamb(lamb) * w * sqrt(nf**2 - ns**2)
N = lambda beta, lamb: beta / k_lamb(lamb) # evaluate N for a given beta, lambda
b_beta_lamb = lambda beta, lamb: (N(beta, lamb)**2 - ns**2) / (nf**2 - ns**2) # evaluate b given beta, lambda
N_from_b = lambda b: abs(sqrt(b * (nf**2 - ns**2) + ns**2)) # evaluate N for a given b value

u0 = 4 * pi * 1e-7
e0 = 8.854e-12
charge = 1.602e-19

divisions=500
plot_divisions=5000
print_N=False
print_neff=False

kf = lambda beta: sqrt(k**2 * nf**2 - beta**2)
ys = lambda beta: sqrt(beta**2 - k**2 * ns**2)
yc = lambda beta: sqrt(beta**2 - k**2 * nc**2)
kf_lamb = lambda beta, lamb: sqrt(k_lamb(lamb)**2 * nf**2 - beta**2)
ys_lamb = lambda beta, lamb: sqrt(beta**2 - k_lamb(lamb)**2 * ns**2)
yc_lamb = lambda beta, lamb: sqrt(beta**2 - k_lamb(lamb)**2 * nc**2)
ys_lamb_2d_x = lambda kx, lamb: sqrt(k_lamb(lamb)**2 * (n1**2 - n2**2) - kx**2)
ys_lamb_2d_y = lambda ky, lamb: sqrt(k_lamb(lamb)**2 * (n1**2 - n2**2) - ky**2)

phi = lambda beta, lamb: atan(yc_lamb(beta, lamb) / kf_lamb(beta, lamb)) # phase constant

default_iter=10

"""
------------------------------------------------------------------------------------------------------------
SOLVERS
------------------------------------------------------------------------------------------------------------
"""
# f_oddNumbers: left hand side         f_evenNumbers: right hand side

# (DEFAULT)
def f1(x):
    return w * kf(x)

def f1_lamb(x, lamb):
    return w * kf_lamb(x, lamb)

def f1_lamb_alt(x, lamb):
    return h * kf_lamb(x, lamb)


def f2(x):
    kappa_f = kf(x)
    return atan(yc(x) / kappa_f) + atan(ys(x) / kappa_f)

def f2_mode(x, m):
    kappa_f = kf(x)
    return atan(yc(x) / kappa_f) + atan(ys(x) / kappa_f) + pi * m

def f2_lamb(x, lamb):
    kappa_f = kf_lamb(x, lamb)
    return atan(yc_lamb(x, lamb) / kappa_f) + atan(ys_lamb(x, lamb) / kappa_f)

def f2_lamb_mode(x, lamb, mode):
    kappa_f = kf_lamb(x, lamb)
    return atan(yc_lamb(x, lamb) / kappa_f) + atan(ys_lamb(x, lamb) / kappa_f) + pi * mode


def solve_default(lamb=wavelength, mode=0, print_N=print_N, max_iter=default_iter):
    """Solves for Beta (TE Mode) under default start and end conditions

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number 0 -> inf
    print_N : bool
        whether or not to print the calculated effective index
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Calculates as a slab waveugide along the x axis (uniform along the y axis)

    """
    start, end = k_lamb(lamb) * ns + 10, k_lamb(lamb) * nf - 10
    x = make_range(start, end, divisions)
    left = [f1_lamb(i, lamb) for i in x]
    right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_lamb(i, lamb) for i in x]
        right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)
    beta = x[stored_index]
    if print_N:
        print(beta / k_lamb(lamb))
    return beta

def solve_default_alt(lamb=wavelength, mode=0, print_N=print_N, max_iter=default_iter):
    """Solves for Beta (TE Mode) under default start and end conditions

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    print_N : bool
        whether or not to print the calculated effective index
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Calculates as a slab waveugide along the y axis (uniform along the x axis)

    """
    start, end = k_lamb(lamb) * ns + 10, k_lamb(lamb) * nf - 10
    x = make_range(start, end, divisions)
    left = [f1_lamb_alt(i, lamb) for i in x]
    right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_lamb_alt(i, lamb) for i in x]
        right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)
    beta = x[stored_index]
    if print_N:
        print(beta / k_lamb(lamb))
    return beta

def solve_range(start, end, lamb=wavelength, mode=0, print_N=print_N, max_iter=default_iter):
    """Solves for Beta (TE Mode) using default solver within some specified start and end range

    PARAMETERS
    ----------
    start : number
        first value to evaluate
    end : number
        last value to evaluate
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    print_N : bool
        whether or not to print the calculated effective index
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Calculates as a slab waveugide along the x axis (uniform along the y axis)

    """
    x = make_range(start, end, divisions)
    left = [f1_lamb(i, lamb) for i in x]
    right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_lamb(i, lamb) for i in x]
        right = [f2_lamb_mode(i, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)
    beta = x[stored_index]
    if print_N:
        print(beta / k_lamb(lamb))
    return beta

# (SOLVING FOR TM MODE)
f1_mag = f1_lamb

def f2_mag(x, lamb, mode):
    kappa_f = kf_lamb(x, lamb)
    return atan((nf / nc)**2 * yc_lamb(x, lamb) / kappa_f) + atan((nf / ns)**2 * ys_lamb(x, lamb) / kappa_f) + pi * mode

def solve_mag(lamb=wavelength, mode=0, print_N=print_N, max_iter=default_iter):
    """Solves for Beta (TM Mode)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    print_N : bool
        whether or not to print the calculated effective index
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Calculates as a slab waveugide along the x axis (uniform along the x axis)

    """
    start, end = k_lamb(lamb) * ns + 10, k_lamb(lamb) * nf - 10
    x = make_range(start, end, divisions)
    left = [f1_mag(i, lamb) for i in x]
    right = [f2_mag(i, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_mag(i, lamb) for i in x]
        right = [f2_mag(i, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)
    beta = x[stored_index]
    if print_N:
        print(beta / k_lamb(lamb))
    return beta

# (EFFECTIVE INDEX METHOD)
def f1_eff(neff, lamb):
    return w * k_lamb(lamb) * sqrt(n1**2 - neff**2)

def f1_eff_alt(neff, lamb):
    return h * k_lamb(lamb) * sqrt(n1**2 - neff**2)

def f3_eff(N, neff, lamb):
    return h * k_lamb(lamb) * sqrt(neff**2 - N**2)

def f3_eff_alt(N, neff, lamb):
    return w * k_lamb(lamb) * sqrt(neff**2 - N**2)


def f2_eff(neff, lamb, mode):
    return mode * pi + atan(sqrt(neff**2 - n3**2) / sqrt(n1**2 - neff**2)) + atan(sqrt(neff**2 - n5**2) / sqrt(n1**2 - neff**2))

def f2_eff_alt(neff, lamb, mode):
    return mode * pi + atan(sqrt(neff**2 - n2**2) / sqrt(n1**2 - neff**2)) + atan(sqrt(neff**2 - n4**2) / sqrt(n1**2 - neff**2))

def f4_eff(N, neff, lamb, mode):
    return mode * pi + atan(sqrt(N**2 - n2**2) / sqrt(neff**2 - N**2)) + atan(sqrt(N**2 - n4**2) / sqrt(neff**2 - N**2))

def f4_eff_alt(N, neff, lamb, mode):
    return mode * pi + atan(sqrt(N**2 - n3**2) / sqrt(neff**2 - N**2)) + atan(sqrt(N**2 - n5**2) / sqrt(neff**2 - N**2))


def solve_eff(lamb=wavelength, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_neff=print_neff, print_N=print_N, max_iter=default_iter):
    """Solves for Beta using the effective index method (considering a rectangular waveguide)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    ret_N : bool
        whether or not to return the calculated N value
    ret_trans : bool
        whether or not to return calculated components kx, ky such that k*nf=√kx^2 + ky^2 + Beta^2
    print_neff : bool
        whether or not to print the calculated neff value
    print_N : bool
        whether or not to print the calculated N value
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number (if !ret_N and !ret_trans)
    
    N : number (if ret_N)
        Calculated effective index

    ky : number (if ret_trans)
        Propagation in y direction

    kx : number (if ret_trans)
        Propagation in x direction

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Assumes a strip waveguide
    • First solves in X direction, then Y

    """
    start, end = max(n2 + 1e-3, n4 + 1e-3), n1 - 1e-3
    x = make_range(start, end, divisions)
    left = [f1_eff(i, lamb) for i in x]
    right = [f2_eff(i, lamb, mode=mode_2) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_eff(i, lamb) for i in x]
        right = [f2_eff(i, lamb, mode=mode_2) for i in x]
        prev_index = find_intersection(left, right)
    neff = x[stored_index]
    neff_stored=neff
    if print_neff:
        print("neff:", neff * k_lamb(lamb))

    start, end = max(n3 + 1e-3, n5 + 1e-3), neff - 1e-3
    x = make_range(start, end, divisions)
    left = [f3_eff(i, neff, lamb) for i in x]
    right = [f4_eff(i, neff, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f3_eff(i, neff, lamb) for i in x]
        right = [f4_eff(i, neff, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)

    if print_N:
        print(x[stored_index]) # x[stored_index] : effective index (N)
    if ret_N and ret_trans:
        return x[stored_index], k_lamb(lamb) * neff_stored, k_lamb(lamb) * x[stored_index] # N, kx, ky
    if ret_N:
        return x[stored_index] # N
    elif ret_trans:
        return k_lamb(lamb) * neff_stored, k_lamb(lamb) * x[stored_index] # kx, ky
    return k_lamb(lamb) * x[stored_index]

def solve_eff_alt(lamb=wavelength, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_neff=print_neff, print_N=print_N, max_iter=default_iter):
    """Solves for Beta using the effective index method (considering a rectangular waveguide)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    ret_N : bool
        whether or not to return the calculated N value
    ret_trans : bool
        whether or not to return calculated components kx, ky such that k*nf=√kx^2 + ky^2 + Beta^2
    print_neff : bool
        whether or not to print the calculated neff value
    print_N : bool
        whether or not to print the calculated N value
    max_iter : number
        maximum number of iterations

    RETURNS
    -------
    Beta : number if (!ret_N and !ret_trans)
    
    N : number (if ret_N)
        Calculated effective index

    kx : number (if ret_trans)
        Propagation in x direction

    ky : number (if ret_trans)
        Propagation in y direction

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Assumes a strip waveguide
    • First solves in Y direction, then X

    """
    start, end = max(n3 + 1e-3, n5 + 1e-3), n1 - 1e-3
    x = make_range(start, end, divisions)
    left = [f1_eff_alt(i, lamb) for i in x]
    right = [f2_eff_alt(i, lamb, mode=mode_2) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f1_eff_alt(i, lamb) for i in x]
        right = [f2_eff_alt(i, lamb, mode=mode_2) for i in x]
        prev_index = find_intersection(left, right)
    neff = x[stored_index]
    if print_neff:
        print("neff:", neff * k_lamb(lamb))
    neff_stored = neff

    start, end = max(n2 + 1e-3, n4 + 1e-3), neff - 1e-3
    x = make_range(start, end, divisions)
    left = [f3_eff_alt(i, neff, lamb) for i in x]
    right = [f4_eff_alt(i, neff, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f3_eff_alt(i, neff, lamb) for i in x]
        right = [f4_eff_alt(i, neff, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)

    if print_N:
        print(x[stored_index]) # x[stored_index] : effective index (N)
    if ret_N and ret_trans:
        return x[stored_index], k_lamb(lamb) * x[stored_index], k_lamb(lamb) * neff_stored # N, kx, ky
    if ret_N:
        return x[stored_index] # N
    if ret_trans:
        return k_lamb(lamb) * x[stored_index], k_lamb(lamb) * neff_stored # kx, ky
    return k_lamb(lamb) * x[stored_index]


# (MARCATILI METHOD)

def f5(kx):
    return w * kx

def f7(ky):
    return h * ky


def f6(kx, lamb, mode):
    return mode * pi + atan(sqrt(k_lamb(lamb)**2 * (n1**2 - n2**2) - kx**2) / kx) + atan(sqrt(k_lamb(lamb)**2 * (n1**2 - n4**2) - kx**2) / kx)

def f8(ky, lamb, mode):
    return mode * pi + atan(sqrt(k_lamb(lamb)**2 * (n1**2 - n3**2) - ky**2) / ky) + atan(sqrt(k_lamb(lamb)**2 * (n1**2 - n5**2) - ky**2) / ky)


def solve_marcatili(lamb=wavelength, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_N=print_N, print_trans=False, max_iter=default_iter, ignore_beta=False):
    """Solves for Beta using the Marcatili method (considering a rectangular waveguide)
       
    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    ret_N : bool
        whether or not to return the calculated N value
    ret_trans : bool
        whether or not to return calculated components kx, ky such that k*nf=√kx^2 + ky^2 + Beta^2
    print_N : bool
        whether or not to print the calculated N value
    print_trans : bool
        whether or not to print kx, ky
    max_iter : number
        maximum number of iterations to run the solver
    ignore_beta : bool
        whether or not to ignore the fact that beta may be less than 0 according to this method

    RETURNS
    -------
    Beta : number if (!ret_N and !ret_trans)
    
    N : number (if ret_N)
        Calculated effective index

    kx : number (if ret_trans)
        Propagation in x direction

    ky : number (if ret_trans)
        Propagation in y direction

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Assumes a strip waveguide

    """

    f_right_1 = f6
    f_right_2 = f8
    start, end = 10, sqrt(k_lamb(lamb)**2 * (n1**2 - ns**2)) - 10
    x = make_range(start, end, divisions)
    left = [f5(i) for i in x]
    right = [f_right_1(i, lamb, mode=mode) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f5(i) for i in x]
        right = [f_right_1(i, lamb, mode=mode) for i in x]
        prev_index = find_intersection(left, right)
    kx = x[stored_index]

    start, end = 10, sqrt(k_lamb(lamb)**2 * (n1**2 - ns**2)) - 10
    x = make_range(start, end, divisions)
    left = [f7(i) for i in x]
    right = [f_right_2(i, lamb, mode=mode_2) for i in x]
    prev_index = find_intersection(left, right)
    if prev_index == None:
        return None
    stored_index = prev_index
    count = 0
    while prev_index != None and count < max_iter:
        count += 1
        start = x[max(prev_index - 5, 0)]
        end = x[min(prev_index + 5, divisions - 1)]
        stored_index = prev_index
        x = make_range(start, end, divisions)
        left = [f7(i) for i in x]
        right = [f_right_2(i, lamb, mode=mode_2) for i in x]
        prev_index = find_intersection(left, right)
    ky = x[stored_index]

    res = (k_lamb(lamb) * n1)**2 - kx**2 - ky**2
    if res < 0 and not ignore_beta:
        return None
    elif res < 0 and ignore_beta:
        return 0
    beta = sqrt(res)
    if print_trans:
        print("kx:", kx, "\tky:", ky)
    if print_N:
        print(beta / k_lamb(lamb))
    if ret_N and ret_trans:
        return beta / k_lamb(lamb), kx, ky # N, kx, ky
    if ret_N:
        return beta / k_lamb(lamb) # N
    if ret_trans:
        return beta, kx, ky
    return beta

def solve_2d_default(lamb=wavelength, mode=0, mode_2=0, ret_trans=False, max_iter=default_iter, print_statements=False):
    """Modified solver of choice to solve for Beta

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    ret_trans : bool
        whether or not to return calculated components kx, ky such that k*nf=√kx^2 + ky^2 + Beta^2
    max_iter : number
        maximum number of iterations to run the solver
    print_statements : bool
        whether or not to print at certain steps of loss calculation

    RETURNS
    -------
    Beta : number
        Propagation in the z direction (propagation constant)
    kx : number
        Propagation in x direction
    ky : number
        Propagation in y direction

    NOTES
    -----
    • If no solution is found (Beta cannot be solved for some mode), it returns None
    • Assumes a strip waveguide

    """
    res = solve_marcatili(lamb=lamb, mode=mode, mode_2=mode_2, ret_trans=True, ignore_beta=True, max_iter=max_iter, print_trans=print_statements, print_N=print_statements)
    if res == None:
        return None
    beta_prev, kx, ky = res

    beta_1 = solve_eff(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_iter, print_neff=print_statements, print_N=print_statements)
    if beta_1 == None:
        return
    beta_2 = solve_eff_alt(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_iter, print_neff=print_statements, print_N=print_statements)
    if beta_2 == None:
        return
    beta = 0.5 * (beta_1 + beta_2)

    temp = sqrt(((k_lamb(lamb) * nf)**2 - (beta)**2) / (kx**2 + ky**2))
    kx *= temp
    ky *= temp

    if ret_trans:
        return beta, kx, ky
    return beta


"""
------------------------------------------------------------------------------------------------------------
PLOTS
------------------------------------------------------------------------------------------------------------
"""

def bV(start=1e-8, end=1e-5, mode=0, iterations=200, solver=solve_default, max_solver_iter=default_iter):
    """Generates values to plot the bV curve

    PARAMETERS
    ----------
    start : number
        first value to evaluate
    end : number
        last value to evaluate
    mode : int
        mode number
    iterations : int
        number of points to evaluate
    solver: func
        solver to use to evaluate beta values
    max_solver_iter : number
        max number of iterations to run the solver

    RETURNS
    -------
    V : list of numbers
        list of V values
    b : list of numbers
        list of b values

    NOTES
    -----
    • Assumes a slab waveguide

    """
    start_lamb = start
    end_lamb = end
    lamb_lst = make_range(start_lamb, end_lamb, iterations, log_opt=True)
    beta_lst = []
    for i in range(len(lamb_lst)):
        l = lamb_lst[i]
        res = solver(lamb=l, mode=mode, max_iter=max_solver_iter)
        beta_lst += [res]
    lamb_lst, beta_lst = remove_none(lamb_lst, beta_lst)
    V = [k_lamb(l) * w * sqrt(nf**2 - ns**2) for l in lamb_lst]
    N_lst = [N(beta_lst[i], lamb_lst[i]) for i in range(len(lamb_lst))]
    b = [(N_lst[i]**2 - ns**2) / (nf**2 - ns**2) for i in range(len(lamb_lst))]
    return V, b

def bV_2d_lamb(start=1e-8, end=1e-5, mode=0, mode_2=0, iterations=200, solver=solve_2d_default, max_solver_iter=default_iter):
    """Generates values to plot the bV curve

    PARAMETERS
    ----------
    start : number
        first value to evaluate
    end : number
        last value to evaluate
    mode : int
        mode number
    iterations : int
        number of points to evaluate
    solver: func
        solver to use to evaluate beta values
    max_solver_iter : number
        max number of iterations to run the solver

    RETURNS
    -------
    V : list of numbers
        list of V values
    b : list of numbers
        list of b values

    NOTES
    -----
    • Assumes a strip waveguide

    """
    start_lamb = start
    end_lamb = end
    lamb_lst = make_range(start_lamb, end_lamb, iterations, log_opt=True)
    beta_lst = []
    for i in range(len(lamb_lst)):
        l = lamb_lst[i]
        res = solver(lamb=l, mode=mode, mode_2=mode_2, max_iter=max_solver_iter)
        beta_lst += [res]
    lamb_lst, beta_lst = remove_none(lamb_lst, beta_lst)
    V = [k_lamb(l) * w * sqrt(nf**2 - ns**2) for l in lamb_lst]
    N_lst = [N(beta_lst[i], lamb_lst[i]) for i in range(len(lamb_lst))]
    b = [(N_lst[i]**2 - ns**2) / (nf**2 - ns**2) for i in range(len(lamb_lst))]
    return V, b

def ey_default(lamb=wavelength, mode=0, margin=1.5, max_solver_iter=default_iter, divisions=plot_divisions, correct_2d=False):
    """Generates values to plot the ey curve along the x-axis given as e_y(x)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        x_vals will be within the range     margin_value * w ± w
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    correct_2d : bool
        whether or not to produce a more accurate field distribution (though some boundary conditions may not be met)

    RETURNS
    -------
    x_vals : list of numbers
        list of x values
    field_vals : list of numbers
        list of e_y(x) computations
    beta : number
        the beta value of the waveguide for a given mode

    NOTES
    -----
    • If no beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Regions (plotting from     -inf -> inf)
        - substrate (left)
        - film (middle)
        - cover (right)
    • List of x_vals returned are in terms of the width

    """
    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * w * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_default(lamb=lamb, mode=mode, max_iter=max_solver_iter)
    if beta == None:
        return None
    if correct_2d:
        s1 = solve_eff(lamb=lamb, mode=mode)
        if s1 == None:
            return
        s2 = solve_eff_alt(lamb=lamb, mode=mode)
        if s2 == None:
            return
        beta_new = 0.5 * (s1+s2)
        res = solve_marcatili(lamb=lamb, mode=mode, ret_trans=True)
        if res == None:
            return
        beta_marc, kx, ky = res
        f_k = sqrt(((k_lamb(lamb) * nf)**2 - (beta_new)**2) / (kx**2 + ky**2))
        factor = sqrt((k_lamb(lamb)*nf)**2 - (f_k * kx)**2) / beta
    else:
        factor = 1
    N = beta / k
    b = (N**2 - ns**2) / (nf**2 - ns**2)

    # -w / 2 center
    start = start_margin - 0.5
    end = end_margin - 0.5
    x_vals = make_range(start, end, divisions)
    i = 0
    field_vals = []
    
    inner = V * sqrt(1 - b)
    back_coeff = sqrt((a + b) / (1 - b))
    
    while x_vals[i] <= -1:
        x = x_vals[i] * w
        exp = V * sqrt(b) * (1 + (x / w))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        field_vals += [val / factor]
        i += 1
        
    while x_vals[i] <= 0:
        x = x_vals[i] * w
        top = V * sqrt(1 - b) * x
        bot = w
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        field_vals += [val * factor]
        i += 1
        
    while i < divisions:
        x = x_vals[i] * w
        top = -V * sqrt(a + b) * x
        bot = w
        exp = top / bot
        val = e**exp
        field_vals += [val / factor]
        i += 1
        
    return x_vals, field_vals, beta

def ey_default_alt(lamb=wavelength, mode=0, margin=1.5, max_solver_iter=default_iter, divisions=plot_divisions, correct_2d=False):
    """Generates values to plot the ey curve along the y-axis given as e_y(y)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    correct_2d : number
        produce a more accurate field distribution (though some boundary conditions may not be met)

    RETURNS
    -------
    y_vals : list of numbers
        list of y values
    field_vals : list of numbers
        list of e_y(y) computations
    beta : number
        the beta value of the waveguide for a given mode

    NOTES
    -----
    • If no beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Regions (plotting from     -inf -> inf)
        - substrate (bottom)
        - film (middle)
        - cover (top)
    • List of y_vals returned are in terms of h

    """
    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * h * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_default_alt(lamb=lamb, mode=mode, max_iter=max_solver_iter)
    if beta == None:
        return None
    if correct_2d:
        s1 = solve_eff(lamb=lamb, mode=mode)
        if s1 == None:
            return
        s2 = solve_eff_alt(lamb=lamb, mode=mode)
        if s2 == None:
            return
        beta_new = 0.5 * (s1+s2)
        res = solve_marcatili(lamb=lamb, mode=mode, ret_trans=True)
        if res == None:
            return
        beta_marc, kx, ky = res
        f_k = sqrt(((k_lamb(lamb) * nf)**2 - (beta_new)**2) / (kx**2 + ky**2))
        factor = sqrt((k_lamb(lamb)*nf)**2 - (f_k * ky)**2) / beta
    else:
        factor = 1
    N = beta / k
    b = (N**2 - ns**2) / (nf**2 - ns**2)

    # -h / 2 center
    start = start_margin - 0.5
    end = end_margin - 0.5
    y_vals = make_range(start, end, divisions)
    i = 0
    field_vals = []
    
    inner = V * sqrt(1 - b)
    back_coeff = sqrt((a + b) / (1 - b))
    
    while y_vals[i] <= -1:
        y = y_vals[i] * h
        exp = V * sqrt(b) * (1 + (y / h))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        field_vals += [val / factor]
        i += 1
        
    while y_vals[i] <= 0:
        y = y_vals[i] * h
        top = V * sqrt(1 - b) * y
        bot = h
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        field_vals += [val * factor]
        i += 1
        
    while i < divisions:
        y = y_vals[i] * h
        top = -V * sqrt(a + b) * y
        bot = h
        exp = top / bot
        val = e**exp
        field_vals += [val / factor]
        i += 1
        
    return y_vals, field_vals, beta

def plot_ey_default(lamb=wavelength, mode=0, margin=1.5, peak_normalize=False, true_normalize=False, center=False, max_solver_iter=default_iter, divisions=plot_divisions, target_file=None, no_output=False):
    """Evaluates the ey values along the x-axis and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        x_vals will be within the range     margin_value * w ± w
    peak_normalize : bool
        whether or not to normalize the plot (adjust Ec such that the plot peaks at 1)
    true_normalize : bool
        whether or not to normalize the plot (adjust Ec such that the plot has total area 1)
    center : bool
        whether or not to center the plot(moving margins of film to -0.5h, 0.5h)
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot

    RETURNS
    -------
    (computed values from ey)
    x_vals : list of numbers
        list of x values
    field_vals : list of numbers
        list of e_y(x) computations
    beta : number
        the beta value of the waveguide for a given mode

    NOTES
    -----
    • If no beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Regions (plotting from     -inf -> inf)
        - substrate (left)
        - film (middle)
        - cover (right)
    • List of x_vals returned are in terms of w

    """
    res = ey_default(lamb=lamb, mode=mode, margin=margin)
    if res == None:
        return # return None
    x_vals, field_vals, beta = res
    
    if peak_normalize == True:
        temp = [abs(val) for val in field_vals]
        temp = max(field_vals)
        ec = 0.99 / temp
        field_vals = [val * ec for val in field_vals]
    elif true_normalize == True:
        a = curve_area(x_vals, field_vals)
        ec = 0.9999/a
        field_vals = [val * ec for val in field_vals]
    i = ((int)(2 * max(field_vals)) + 1) / 2
    if center == True:
        x_vals = make_range(-margin, margin, divisions)
        plt.axvline(-0.5, color="red", linestyle="--")
        plt.axvline(0.5, color="red", linestyle="--")
        if min(field_vals) < -0.1:
            plt.axis([-margin, margin, -i, i])
        else:
            plt.axis([-margin, margin, 0, i])
    else:
        plt.axvline(-1, color="red", linestyle="--")
        plt.axvline(0, color="red", linestyle="--")
        if min(field_vals) < -0.5:
            plt.axis([-margin - 0.5, margin - 0.5, -i, i])
        else:
            plt.axis([-margin - 0.5, margin - 0.5, 0, i])
        
    plt.plot(x_vals, field_vals)

    if target_file != None and type(target_file) == str:
        savefig(target_file)

    return x_vals, field_vals, beta

def ey_2d_eff(lamb=wavelength, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=default_iter, divisions=plot_divisions):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    margin_x : number
        x_vals will be within the range     margin_value * w ± w
    margin_y : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate

    RETURNS
    -------
    (computed values from ey)
    x_vals : list of numbers
        list of x values
    x_field_vals : list of numbers
        list of e_y(x) computations
    y_vals : list of numbers
        list of y values
    y_field_vals : list of numbers
        list of e_y(y) computations
    beta : number
        the Beta value of the waveguide for a given mode

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide
    • Uses the effective index method
    • List of x_vals, y_vals returned are in terms of w, h

    """
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    res = solve_eff(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_solver_iter, ret_trans=True)
    if res == None:
        return None
    beta_prev, beta=res
    N = beta / k
    N_prev = beta_prev / k

    start_margin=-margin_y
    end_margin=margin_y
    start = start_margin - 0.5
    end = end_margin - 0.5

    # Y REGION
    y_vals = make_range(start, end, divisions)
    V = k * h * sqrt(N_prev**2 - n5**2)
    a = (n5**2 - n3**2) / (N_prev**2 - n5**2)
    b = (N**2 - n5**2) / (N_prev**2 - n5**2)

    i = 0
    y_field_vals = []
    inner = V * sqrt(1 - b)
    back_coeff = sqrt((a + b) / (1 - b))
    
    while y_vals[i] <= -1:
        y = y_vals[i] * h
        exp = V * sqrt(b) * (1 + (y / h))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        y_field_vals += [val]
        i += 1
        
    while y_vals[i] <= 0:
        y = y_vals[i] * h
        top = V * sqrt(1 - b) * y
        bot = h
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        y_field_vals += [val]
        i += 1
        
    while i < divisions:
        y = y_vals[i] * h
        top = -V * sqrt(a + b) * y
        bot = h
        exp = top / bot
        val = e**exp
        y_field_vals += [val]
        i += 1

    # X REGION
    start_margin=-margin_x
    end_margin=margin_x
    start = start_margin - 0.5
    end = end_margin - 0.5

    x_vals = make_range(start, end, divisions)
    V = k * w * sqrt(n1**2 - n4**2)
    a = (n4**2 - n2**2) / (n1**2 - n4**2)
    b = (N_prev**2 - n4**2) / (n1**2 - n4**2)
    d = (n2 / n1)**2

    i = 0
    x_field_vals = []
    inner = V * sqrt(1 - b)
    back_coeff = (1 / d) * sqrt((a + b) / (1 - b))
    
    while x_vals[i] <= -1:
        x = x_vals[i] * w
        exp = V * sqrt(b) * (1 + (x / w))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        x_field_vals += [val]
        i += 1
        
    while x_vals[i] <= 0:
        x = x_vals[i] * w
        top = V * sqrt(1 - b) * x
        bot = w
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        x_field_vals += [val]
        i += 1
        
    while i < divisions:
        x = x_vals[i] * w
        top = -V * sqrt(a + b) * x
        bot = w
        exp = top / bot
        val = e**exp
        x_field_vals += [val]
        i += 1

    return x_vals, x_field_vals, y_vals, y_field_vals, beta

def ey_2d_eff_alt(lamb=wavelength, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=default_iter, divisions=plot_divisions):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    margin_x : number
        x_vals will be within the range     margin_value * w ± w
    margin_y : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate

    RETURNS
    -------
    (computed values from ey)
    x_vals : list of numbers
        list of x values
    x_field_vals : list of numbers
        list of e_y(x) computations
    y_vals : list of numbers
        list of y values
    y_field_vals : list of numbers
        list of e_y(y) computations
    beta : number
        the Beta value of the waveguide for a given mode

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide
    • Uses the alternate effective index method
    • List of x_vals, y_vals returned are in terms of w, h

    """
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    res = solve_eff_alt(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_solver_iter, ret_trans=True)
    if res == None:
        return None
    beta, beta_prev=res
    N = beta / k
    N_prev = beta_prev / k

    start_margin=-margin_x
    end_margin=margin_x
    start = start_margin - 0.5
    end = end_margin - 0.5

    # Y REGION
    y_vals = make_range(start, end, divisions)
    V = k * h * sqrt(n1**2 - n5**2)
    a = (n5**2 - n3**2) / (n1**2 - n5**2)
    b = (N_prev**2 - n5**2) / (n1**2 - n5**2)
    # print(V, a, b)

    i = 0
    y_field_vals = []
    inner = V * sqrt(1 - b)
    back_coeff = sqrt((a + b) / (1 - b))
    
    while y_vals[i] <= -1:
        y = y_vals[i] * h
        exp = V * sqrt(b) * (1 + (y / h))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        y_field_vals += [val]
        i += 1
        
    while y_vals[i] <= 0:
        y = y_vals[i] * h
        top = V * sqrt(1 - b) * y
        bot = h
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        x_field_vals += [val]
        i += 1
        
    while i < divisions:
        y = y_vals[i] * h
        top = -V * sqrt(a + b) * y
        bot = h
        exp = top / bot
        val = e**exp
        y_field_vals += [val]
        i += 1

    # X REGION
    start_margin=-margin_x
    end_margin=margin_x
    start = start_margin - 0.5
    end = end_margin - 0.5

    x_vals = make_range(start, end, divisions)
    V = k * w * sqrt(N_prev**2 - n4**2)
    a = (n4**2 - n2**2) / (N_prev**2 - n4**2)
    b = (N**2 - n4**2) / (N_prev**2 - n4**2)
    d = (n2 / N_prev)**2

    i = 0
    x_field_vals = []
    inner = V * sqrt(1 - b)
    back_coeff = (1 / d) * sqrt((a + b) / (1 - b))
    
    while x_vals[i] <= -1:
        x = x_vals[i] * w
        exp = V * sqrt(b) * (1 + (x / w))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        y_field_vals += [val]
        i += 1
        
    while x_vals[i] <= 0:
        x = x_vals[i] * w
        top = V * sqrt(1 - b) * x
        bot = w
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        x_field_vals += [val]
        i += 1
        
    while i < divisions:
        x = x_vals[i] * w
        top = -V * sqrt(a + b) * x
        bot = w
        exp = top / bot
        val = e**exp
        x_field_vals += [val]
        i += 1

    return x_vals, x_field_vals, y_vals, y_field_vals, beta

def ey_2d_marcatili(lamb=wavelength, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=default_iter, divisions=plot_divisions, ignore_beta=False):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    margin_x : number
        x_vals will be within the range     margin_value * w ± w
    margin_y : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate

    RETURNS
    -------
    (computed values from ey)
    x_vals : list of numbers
        list of x values
    x_field_vals : list of numbers
        list of e_y(x) computations
    y_vals : list of numbers
        list of y values
    y_field_vals : list of numbers
        list of e_y(y) computations
    beta : number
        the Beta value of the waveguide for a given mode

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide
    • Uses the marcatili method
    • List of x_vals, y_vals returned are in terms of w, h

    """
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    res = solve_marcatili(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_solver_iter, ret_trans=True, ignore_beta=ignore_beta)
    if res == None:
        return None
    beta_prev, kx, ky = res
    beta = beta_prev

    start_margin=-margin_x
    end_margin=margin_x
    start = start_margin - 0.5
    end = end_margin - 0.5

    beta_x = sqrt((k_lamb(lamb)*n1)**2 - kx**2)
    beta_y = sqrt((k_lamb(lamb)*n1)**2 - ky**2)
    
    # X REGION
    x_vals = make_range(start, end, divisions)
    yc_x = ys_x = ys_lamb_2d_x(kx, lamb)

    i = 0
    x_field_vals = []
    constant = cos(kx * w) + (yc_x / kx) * sin(kx * w)
    
    while x_vals[i] <= -1:
        x = x_vals[i] * w
        val = constant * e**(ys_x * (x + w))
        x_field_vals += [val]
        i += 1

    while x_vals[i] <= 0:
        x = x_vals[i] * w
        val = cos(kx * x) - (yc_x / kx) * sin(kx * x)
        x_field_vals += [val]
        i += 1

    while i < divisions:
        x = x_vals[i] * w
        val = e**(-x * yc_x)
        x_field_vals += [val]
        i += 1

    # Y REGION
    start_margin=-margin_y
    end_margin=margin_y
    start = start_margin - 0.5
    end = end_margin - 0.5

    y_vals = make_range(start, end, divisions)
    yc_y = ys_y = ys_lamb_2d_y(ky, lamb)
    constant = cos(ky * h) + (yc_y / ky) * sin(ky * h)

    i = 0
    y_field_vals = []

    while y_vals[i] <= -1:
        y = y_vals[i] * h
        val = constant * e**(ys_y * (y + h))
        y_field_vals += [val]
        i += 1
        
    while y_vals[i] <= 0:
        y = y_vals[i] * h
        val = cos(ky * y) - (yc_y / ky) * sin(ky * y)
        y_field_vals += [val]
        i += 1
        
    while i < divisions:
        y = y_vals[i] * h
        val = e**(-y * yc_y)
        y_field_vals += [val]
        i += 1

    return x_vals, x_field_vals, y_vals, y_field_vals, beta

def plot_2d(ey=ey_2d_marcatili, lamb=wavelength, mode=0, mode_2=0, margin_x=2, margin_y=2, max_solver_iter=default_iter, divisions=plot_divisions, normalize=False, target_file=None):
    """Plots the mode profile looking at a cross-section of the waveguide

    PARAMETERS
    ----------
    ey : func
        chosen function to evaluate ey values
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    margin_x : number
        x_vals will be within the range     margin_value * w ± w
    margin_y : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    normalize : bool
        whether or not to normalize the power to 1W
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot

    RETURNS
    -------
    (computed values from ey)
    x_vals : list of numbers
        list of x values
    x_field_vals : list of numbers
        list of e_y(x) computations
    y_vals : list of numbers
        list of y values
    y_field_vals : list of numbers
        list of e_y(y) computations
    beta : number
        the Beta value of the waveguide for a given mode

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide
    • Uses the marcatili method
    • List of x_vals, y_vals returned are in terms of w, h

    """
    res = ey(lamb=lamb, mode=mode, mode_2=mode_2, margin_x=margin_x, margin_y=margin_y, max_solver_iter=max_solver_iter, divisions=divisions)
    if res == None:
        return
    if normalize:
        a_x = curve_area(res[0], [i**2 for i in res[1]], total_area=True)
        res1 = [i / a_x for i in res[1]]
        a_y = curve_area(res[2], [i**2 for i in res[3]], total_area=True)
        res3 = [i / a_y for i in res[3]]
    else:
        res1 = res[1]
        res3 = res[3]

    x = np.array([res1])
    y = np.array([res3])
    r = np.dot(abs(y.T), abs(x))

    plt.imshow(r).set_extent([res[0][0], res[0][divisions - 1], res[2][0], res[2][divisions - 1]])
    plt.colorbar()
    plt.xticks([-1, 0])
    plt.yticks([-1, 0])
    plt.axvline(-1, color="white", linestyle="--")
    plt.axvline(0, color="white", linestyle="--")
    plt.axhline(-1, color="white", linestyle="--")
    plt.axhline(0, color="white", linestyle="--")
    plt.xlabel("x")
    plt.ylabel("y")

    if target_file != None and type(target_file) == str:
        savefig(target_file)

    return res

def em_default(lamb=wavelength, mode=0, margin=1.5, max_solver_iter=default_iter, divisions=divisions):
    """Generates values to plot the TM mode along the x-axis given as h_y(x)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        x_vals will be within the range     margin_value * w ± w
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate

    RETURNS
    -------
    x_vals : list of numbers
        list of x values
    field_vals : list of numbers
        list of h_y(x) computations
    beta : number
        the beta value of the waveguide for a given mode

    NOTES
    -----
    • If no beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Regions (plotting from     -inf -> inf)
        - substrate (left)
        - film (middle)
        - cover (right)
    • List of x_vals returned are in terms of w

    """
    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * w * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_mag(lamb=lamb, mode=mode, max_iter=max_solver_iter)
    if beta == None:
        return None
    N = beta / k
    b = (N**2 - ns**2) / (nf**2 - ns**2)
    d = (nc / nf)**2

    # -w / 2 center
    start = start_margin - 0.5
    end = end_margin - 0.5
    x_vals = make_range(start, end, divisions)
    i = 0
    field_vals = []
    
    inner = V * sqrt(1 - b)
    back_coeff = sqrt((a + b) / (1 - b)) / d
    
    while x_vals[i] <= -1:
        x = x_vals[i] * w
        exp = V * sqrt(b) * (1 + (x / w))
        val = (cos(inner) + sin(inner) * back_coeff) * e**exp
        field_vals += [val]
        i += 1
    
    while x_vals[i] <= 0:
        x = x_vals[i] * w
        top = V * sqrt(1 - b) * x
        bot = w
        front = cos(top / bot)
        back = sin(top / bot)
        val = front - back * back_coeff
        field_vals += [val]
        i += 1
        
    while i < divisions:
        x = x_vals[i] * w
        top = -V * sqrt(a + b) * x
        bot = w
        exp = top / bot
        val = e**exp
        field_vals += [val]
        i += 1
        
    return x_vals, field_vals, beta

def plot_em_default(target_file=None):
    """Plots the TM mode
    PARAMETERS
    ----------
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot

    RETURNS
    -------
    x_vals : list of numbers
        list of x values
    field_vals : list of numbers
        list of h_y(x) computations
    beta : number
        the beta value of the waveguide for a given mode

    NOTES
    -----
    • If no beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Regions (plotting from     -inf -> inf)
        - substrate (left)
        - film (middle)
        - cover (right)
    • This is the basic plotting code. It is suggested you evaluate and plot em mode yourself using em_default

    """
    res = em_default()
    if res == None:
        return 
    plt.plot(res[0], res[1])
    if target_file != None and type(target_file) == str:
        savefig(target_file)
    return res

"""
------------------------------------------------------------------------------------------------------------
LOSS FUNCTIONS
------------------------------------------------------------------------------------------------------------
"""

# MATH FUNCTIONS
rad = lambda x: x * pi / 180
scalar = lambda x: x * 180 / pi

def absorption_loss(N=1e18, lamb=wavelength, doping="n", print_afc=False, constant=False):
    """Calculates the loss due to free carrier absorption
    PARAMETERS
    ----------
    N : number
        doping concentration
    lamb : number
        wavelength (lambda)
    doping : string
        doping type (either 'n' or 'p')
    print_afc : bool
        whether or not to print per unit dopant concentration
    constant : bool
        whether or not loss does not vary with effective index contrast

    RETURNS
    -------
    loss : number
        Absorption loss due to free carrier (intraband) absorption

    NOTES
    -----
    • Assumes that the material waveguide used is heavily doped silicon
    • If it is neither 'n' nor 'p' type doped, it returns None

    """
    if constant == True:
        return absorption_default(N, lamb, doping, print_afc)
    if doping == "n":
        m = 0.26 * 9.10938356e-31
        u = 1318 / (1 + (N / 1e17)**0.85) + 92
    elif doping == "p":
        m = 0.39 * 9.10938356e-31
        u = 420 / (1 + (N / 1.6e17)**0.7) + 50
    else:
        return # ERROR: can't be neither typed

    top = charge**3 * lamb**2 * 1e8
    bot = 4 * pi**2 * nf * m**2 * e0 * u * c_light**3
    afc = top / bot

    if print_afc:
        print(afc)
    return afc * N

def absorption_default(N=1e18, lamb=wavelength, doping="n", print_afc=False):
    """Caluclates the Absorption Loss assuming the film has a refractive index of 3.5

    PARAMETERS
    ----------
    N : integer
        Doping Concentration
    lamb : number
        wavelength (lambda)
    doping : string
        doping type
    print_afc : bool
        whether or not to print per unit dopant concentration

    RETURNS
    -------
    loss : number
        Absorption loss due to free carrier (intraband) absorption

    NOTES
    -----
    • Assumes that the material waveguide used is heavily doped silicon
    • If it is neither 'n' nor 'p' type doped, it returns None

    """
    if doping == "n":
        m = 0.26 * 9.10938356e-31
        u = 1318 / (1 + (N / 1e17)**0.85) + 92
    elif doping == "p":
        m = 0.39 * 9.10938356e-31
        u = 420 / (1 + (N / 1.6e17)**0.7) + 50
    else:
        return # ERROR: can't be neither typed

    top = charge**3 * lamb**2 * 1e8
    bot = 4 * pi**2 * 3.5 * m**2 * e0 * u * c_light**3
    afc = top / bot

    if print_afc:
        print(afc)
    return afc * N

def roughness_loss(lamb=wavelength, mode=0, roughness=11.5e-9, print_statements=False):
    """Calculates the loss due to surface roughness

    Parameters
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    roughness : number
        roughness variance
    print_statements : bool
        whether or not to print at certain steps of loss calculation

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide along the x-axis (uniform along the y-axis)

    """
    res = propagation_angle(lamb=lamb, mode=mode, ret_beta=True)
    if res == None:
        return
    theta_incident, beta = res
    if theta_incident == None:
        return
    theta_m = 90 - theta_incident
    A = 4 * pi * roughness / lamb
    if print_statements:
        print("theta", theta_incident)
        print("A", A)

    return A**2 * (cos(rad(theta_m)))**3 / (2 * sin(rad(theta_m)) * effective_film_thickness_x(beta, lamb)) * 4.343 / 100

def roughness_loss_2d(lamb=wavelength, mode=0, mode_2=0, roughness=11e-9, correct_theta=False, print_statements=False, method=1, fixed_thickness=False):
    """Calculates the loss due to surface roughness assuming a 2D waveguide

    Parameters
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    roughness : number
        roughness variance
    correct_theta : bool
        whether to use a larger propagation angle (treating it similar to a slab waveguide)
    print_statements : bool
        whether or not to print at certain steps of loss calculation
    method : opt
        OPT     METHOD
        0       Directly sum attenuation from interaction with horizontal and vertical sidewalls
        1       Sum the squares of components of loss from interaction with horizontal and vertical sidewalls eg. sqrt(a_x^2 + a_y^2)
        2       Return the seperate components of loss from interaction with horizontal and vertical sidewalls
    fixed_thickness : bool

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide

    """
    res = solve_marcatili(lamb=wavelength, mode=mode, mode_2=mode_2, ret_trans=True)
    if res == None:
        return
    beta_prev, kx, ky = res
    beta_1 = solve_eff(lamb=lamb, mode=mode, mode_2=mode_2)
    if beta_1 == None:
        return
    beta_2 = solve_eff_alt(lamb=lamb, mode=mode, mode_2=mode_2)
    if beta_2 == None:
        return
    beta = 0.5 * (beta_1 + beta_2)
    temp = sqrt(((k_lamb(lamb) * nf)**2 - (beta)**2) / (kx**2 + ky**2))
    kx *= temp
    ky *= temp
    A = 4 * pi * roughness / lamb
    if correct_theta:
        theta_x = scalar(atan(kx/sqrt(beta**2 + ky**2)))
        theta_y = scalar(atan(ky/sqrt(beta**2 + kx**2)))
    else:
        theta_x = scalar(atan(kx/beta))
        theta_y = scalar(atan(ky/beta))
    theta_m_x = 90 - theta_x
    theta_m_y = 90 - theta_y

    if fixed_thickness:
        thickness_x = fixed_w
        thickness_y = fixed_h
    else:
        thickness_x = effective_film_thickness_x(kx, lamb)
        thickness_y = effective_film_thickness_y(ky, lamb)

    loss_x = (cos(rad(theta_m_x)))**3 / (2 * sin(rad(theta_m_x)) * thickness_x)
    loss_y = (cos(rad(theta_m_y)))**3 / (2 * sin(rad(theta_m_y)) * thickness_y)
    if print_statements:
        print("beta:", beta, "kx:", kx, "ky", ky)
        print("theta_x", theta_x, "theta_y", theta_y)
        print("theta_m_x:", theta_m_x, "theta_m_y:", theta_m_y)
        print("thickness_x:", thickness_x, "thickness_y:", thickness_y)
        print("loss_components", A**2 * loss_x * 4.343 / 100, A**2 * loss_y * 4.343 / 100)
        print("per unit roughness", A**2 * (loss_x + loss_y) * 4.343 / 100 / roughness**2)
    if method == 0:
        return A**2 * (loss_x + loss_y) * 4.343 / 100
    elif method == 1:
        return A**2 * sqrt(loss_x**2 + loss_y**2) * 4.343 / 100
    elif method == 2:
        return A**2 * loss_x * 4.343 / 100, A**2 * loss_y * 4.343 / 100

def roughness_loss_solver(lamb=wavelength, mode=0, roughness=11e-9, solver=solve_default, print_statements=False):
    """Calculates the loss due to surface roughness given some solver

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    roughness : number
        roughness variance
    print_statements : bool
        whether or not to print at certain steps of loss calculation

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide
    • Uses the method established by Tien

    """
    beta = solver(lamb=lamb, mode=mode)
    if beta == None:
        return
    theta_incident = scalar(atan(kf_lamb(beta, lamb) / beta))
    theta_m = 90 - theta_incident
    A = 4 * pi * roughness / lamb
    if print_statements:
        print("beta:", beta)
        print("theta_m:", theta_m)
        print("A:", A)
    return A**2 * (cos(rad(theta_m)))**3 / (2 * sin(rad(theta_m)) * effective_film_thickness_x(beta, lamb)) * 4.343 / 100

def roughness_loss_theta(theta_incident=23, lamb=wavelength, mode=0, roughness=11e-9, solver=solve_default, print_statements=False):
    """Calculates the loss due to surface roughness given some angle of propagation

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    roughness : number
        roughness variance
    solver : func
        method to solve for beta
    print_statements : bool
        whether or not to print at certain steps of loss calculation

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide

    """
    theta_m = 90 - theta_incident
    beta = solver(lamb, mode)
    if beta == None:
        return
    A = 4 * pi * roughness / lamb
    if print_statements:
        print("beta:", beta)
        print("theta_m:", theta_m)
        print("A:", A)
    return A**2 * (cos(rad(theta_m)))**3 / (2 * sin(rad(theta_m)) * effective_film_thickness_x(beta, lamb)) * 4.343 / 100

"""
------------------------------------------------------------------------------------------------------------
ANGLES
------------------------------------------------------------------------------------------------------------
"""

def propagation_angle(lamb=wavelength, mode=0, ret_beta=False):
    """Computes the angle of propagation within a slab waveguide

    ----------------------- top cladding
      θ /\        /\
       /  \      /  \
      /    \    /    \    / waveguide core 
     /      \  /      \  /        
    /        \/        \/
    ----------------------- bottom cladding
    -> light propagation
    

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    ret_beta : bool
        whether or not to return beta

    RETURNS
    -------
    theta : number
        propagation angle

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide

    """
    beta = solve_default(lamb, mode)
    if beta == None:
        return None
    theta = scalar(atan(kf_lamb(beta, lamb) / beta))
    if ret_beta:
        return theta, beta
    return theta

def propagation_angle_2d(lamb=wavelength, mode=0, ret_beta=False):
    """Computes the angle of propagation within a strip waveguide

    ----------------------- top cladding
      θ /\        /\
       /  \      /  \
      /    \    /    \    / waveguide core 
     /      \  /      \  /        
    /        \/        \/
    ----------------------- bottom cladding
    -> light propagation
    

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    ret_beta : bool
        whether or not to return beta

    RETURNS
    -------
    theta_x : number
        propagation angle in the x direction
    theta_y : number
        propagation angle in the y direction

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide

    """
    res = solve_marcatili(lamb=lamb, mode=mode, ret_trans=True)
    if res == None:
        return None
    beta, kx, ky = res
    if ret_beta:
        return beta, scalar(atan(kx / beta)), scalar(atan(ky / beta))
    return scalar(atan(kx / beta)), scalar(atan(ky / beta))

def critical_angle():
    """Calculates the critical angle (between light and the vector normal to the waveguide boundary) of a slab waveguide
    ------------ top cladding
        /|\ 
       /θ| \
      /  |  \    waveguide core 
     /   |   \  
    /    |    \
    ------------ bottom cladding
    -> light propagation

    RETURNS
    -------
    theta : number
        critical angle (see picture)

    """
    return scalar(asin(ns / nf))

def critical_angle_incident():
    """Calculates the critical angle (between light and the vector normal to the waveguide boundary) of a slab waveguide
    ----------------------top cladding
      θ /|\ 
       / | \
      /  |  \
     /   |   \  waveguide core 
    /    |    \
    
    -> light propagation
    ----------------------bottom cladding

    RETURNS
    -------
    theta : number
        critical angle (see picture)

    """
    return 90 - critical_angle()

"""
------------------------------------------------------------------------------------------------------------
EFFECTIVE FILM THICKNESS
------------------------------------------------------------------------------------------------------------
"""

def effective_film_thickness_default(mode=0, lamb=1550e-9, solver=solve_default):
    """Computes the effective thickness of the waveguide (width of the film)

    PARAMETERS
    ----------
    beta : number
        
    lamb : number
        wavelength (lambda)
    solver : func
        solver of choice to calculate beta

    RETURNS
    -------
    thickness : number
        effective width of the film

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab waveguide

    """
    beta = solver(mode=mode, lamb=lamb)
    if beta == None:
        return
    return effective_film_thickness_x(beta, lamb)

def effective_film_thickness_x(kx, lamb):
    """Computes the effective width of the waveguide film


    PARAMETERS
    ----------
    kx : number
        
    lamb : number
        wavelength (lambda)
    solver : func
        solver of choice to calculate beta

    RETURNS
    -------
    thickness : number
        effective width of the film

    NOTES
    -----
    • Assumes a strip waveguide

    """
    if kx**2 >= k_lamb(lamb)**2 * (n1**2 - n2**2):
        return w
    return w + 2 / ys_lamb_2d_x(kx, lamb)

def effective_film_thickness_y(ky, lamb):
    """Computes the effective height of the waveguide film

    PARAMETERS
    ----------
    beta : number
        
    lamb : number
        wavelength (lambda)
    solver : func
        solver of choice to calculate beta

    RETURNS
    -------
    thickness : number
        effective height of the film

    NOTES
    -----
    • Assumes a strip waveguide

    """
    if ky**2 >= k_lamb(lamb)**2 * (n1**2 - n2**2):
        return h
    return h + 2 / ys_lamb_2d_y(ky, lamb)

def effective_film_thickness_2d(lamb=wavelength, mode=0, mode_2=0, fixed_thickness=False):
    """computes the effective thickness of the waveguide in the x and y direction

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in y direction
    mode_2 : int
        mode number for propagation in x direction
    fixed_thickness : bool
        whether or not to neglect penetration depth

    RETURNS
    -------
    thickness_x : number
        effective width of the film
    thickness_y : number
        effective height of the film

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide

    """
    if fixed_thickness:
        return w, h

    res = solve_marcatili(lamb=lamb, mode=mode, mode_2=mode_2, ret_trans=True, ignore_beta=True)
    if res == None:
        return None
    beta_prev, kx, ky = res

    beta_1 = solve_eff(lamb=lamb, mode=mode, mode_2=mode_2)
    if beta_1 == None:
        return
    beta_2 = solve_eff_alt(lamb=lamb, mode=mode, mode_2=mode_2)
    if beta_2 == None:
        return
    beta = 0.5 * (beta_1 + beta_2)

    temp = sqrt(((k_lamb(lamb) * nf)**2 - (beta)**2) / (kx**2 + ky**2))
    kx *= temp
    ky *= temp

    thickness_x = effective_film_thickness_x(kx, lamb)
    thickness_y = effective_film_thickness_y(ky, lamb)
    return thickness_x, thickness_y  # thickness_x (w), thickness_y (h)

"""
------------------------------------------------------------------------------------------------------------
CONFINEMENT FACTOR
------------------------------------------------------------------------------------------------------------
"""

def confinement_factor_single_mode(lamb=wavelength, solver=solve_default):
    """Computes the confinement factor looking only at the zero-th mode (m = 0)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    solver : func
        method to solve for beta

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a slab waveguide along the x axis

    """
    beta = solver(lamb=lamb)
    yc_val = yc_lamb(beta, lamb)
    ys_val = ys_lamb(beta, lamb)
    kf_val = kf_lamb(beta, lamb)
    top = w + yc_val / (kf_val**2 + yc_val**2) + ys_val / (kf_val**2 + ys_val**2)
    bot = w + 1 / yc_val + 1 / ys_val
    return top / bot

def confinement_factor(lamb=wavelength, mode=0, margin=20, print_statements=False):
    """Computes the confinement factor looking at a single mode

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        x_vals will be within the range     margin_value * w ± w
    print_statements : bool
        whether or not to print at certain steps of calculation

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a slab waveguide along the x axis

    """
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_default(lamb=lamb, mode=mode, margin=margin)
    total = 0
    film = 0
    if res == None:
        return
    x_vals, field_vals, beta = res
    x_vals = [i * w for i in x_vals]
    coeff = beta / (2 * omega * u0)
    ec = 1
    phase = phi(beta, lamb)
    ef = ec / cos(phase)
    es = ef * cos(phase - w * kf_lamb(beta, lamb))
    f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + w))
    total += integrate(x_vals, end=-w, func=f1)
    f2 = lambda x: coeff * ef**2 * (cos(kf_lamb(beta, lamb) * x + phase))**2
    temp = integrate(x_vals, start=-w, end=0, func=f2)
    film += temp
    total += temp
    f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
    total += integrate(x_vals, start=0, func=f3)
    if print_statements:
        print("film_power", film, "total_power", total)
    return film / total

def confinement_factor_alt(lamb=wavelength, mode=0, margin=20, print_statements=False):
    """Computes the confinement factor looking at a single mode

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        y_vals will be within the range     margin_value * h ± h
    print_statements : bool
        whether or not to print at certain steps of calculation

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a slab waveguide along the y axis

    """
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_default_alt(lamb=lamb, mode=mode, margin=margin)
    total = 0
    film = 0
    if res == None:
        return
    x_vals, field_vals, beta = res
    x_vals = [i * h for i in x_vals]
    coeff = beta / (2 * omega * u0)
    ec = 1
    phase = phi(beta, lamb)
    ef = ec / cos(phase)
    es = ef * cos(phase - h * kf_lamb(beta, lamb))
    f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + h))
    total += integrate(x_vals, end=-h, func=f1)
    f2 = lambda x: coeff * ef**2 * (cos(kf_lamb(beta, lamb) * x + phase))**2
    temp = integrate(x_vals, start=-h, end=0, func=f2)
    film += temp
    total += temp
    f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
    total += integrate(x_vals, start=0, func=f3)
    if print_statements:
        print("film_power", film, "total_power", total)
    return film / total

def confinement_factor_interface(roughness=11.5e-9, lamb=wavelength, mode=0, margin=20, print_statements=False):
    """Computes the confinement factor in the cladding looking at a single mode

    PARAMETERS
    ----------
    roughness : number
        roughness variance denoting the range surrounding the cladding to consider
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        x_vals will be within the range     margin_value * w ± w
    print_statements : bool
        whether or not to print at certain steps of calculation

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a slab waveguide along the x axis

    """
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_default(lamb=lamb, mode=mode, margin=margin)
    total = 0
    film = 0
    if res == None:
        return
    x_vals, field_vals, beta = res
    clad = curve_area(x_vals, field_vals, start = -1 - roughness / w, end = -1 + roughness / w)
    total = curve_area(x_vals, field_vals)
    if print_statements:
        print("interface_power", clad, "total_power", total)
    return 2 * clad / total

def confinement_factor_interface_2d(roughness=11.5e-9, lamb=wavelength, mode=0, mode_2=0, margin=20, divisions=plot_divisions, include_corners=True, correct_2d=True):
    """Computes the confinement factor at the waveguide boundary interface looking at a single mode

    PARAMETERS
    ----------
    roughness : number
        roughness variance denoting the range surrounding the cladding to consider
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    margin : number
        x_vals will be within the range     margin_value * w ± w
    divisions : number
        number of points to evaluate
    correct_2d : bool
        whether or not to produce a more accurate field distribution (though some boundary conditions may not be met)
    include_corners : bool
        whether or not to include the field of the corners in confinement calculation

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a strip waveguide

    """
    divisions = (divisions) * 1 * margin
    freq = c_light / lamb
    omega = freq * 2 * pi
    if mode_2 != 0:
        mod_correct = False
    else:
        mod_correct=correct_2d
    res = ey_default(lamb=lamb, mode=mode, margin=margin, divisions=divisions, correct_2d=mod_correct)
    if res == None:
        return
    x_vals, x_field_vals, beta = res
    if mode_2 != 0:
        mod_correct = False
    else:
        mod_correct=correct_2d
    res = ey_default_alt(lamb=lamb, mode=mode_2, margin=margin, divisions=divisions, correct_2d=mod_correct)
    if res == None:
        return res
    y_vals, y_field_vals, beta = res

    dx = (2 * margin * w) / divisions
    dy = (2 * margin * h) / divisions

    x_start = 0
    while x_vals[x_start] <= -1:
        x_start += 1
    x_end = x_start
    while x_vals[x_end] <= 0:
        x_end += 1
    x_end -= 1
    y_start = 0
    while y_vals[y_start] < -1:
        y_start += 1
    y_end = y_start
    while y_vals[y_end] < 0:
        y_end += 1
    y_end -= 1

    # AIR
    x_cur = x_start
    middle_power_x = 0
    while x_cur < x_end:
        middle_power_x += x_field_vals[x_cur]**2
        x_cur += 1
    middle_power_x *= dx
    y_cur = 0
    n2_power = 0
    while y_cur < y_start:
        n2_power += middle_power_x * y_field_vals[y_cur]**2
        y_cur += 1
    n2_power *= dy

    # LEFT MIDDLE SUBSTRATE
    x_cur = 0
    left_power_x = 0
    while x_cur < x_start:
        left_power_x += x_field_vals[x_cur]**2
        x_cur += 1
    left_power_x *= dx
    n5_power = 0
    while y_cur < y_end:
        n5_power += left_power_x * y_field_vals[y_cur]**2
        y_cur += 1
    n5_power *= dy

    # FILM
    y_cur = y_start
    n1_power = 0
    while y_cur < y_end:
        n1_power += middle_power_x * y_field_vals[y_cur]**2
        y_cur += 1
    n1_power *= dy

    # RIGHT MIDDLE SUBSTRATE
    x_cur = x_end
    right_power_x = 0
    while x_cur < divisions:
        right_power_x += x_field_vals[x_cur]**2
        x_cur += 1
    right_power_x *= dx
    y_cur = y_start
    n3_power = 0
    while y_cur < y_end:
        n3_power += right_power_x * y_field_vals[y_cur]**2
        y_cur += 1
    n3_power *= dy

    # BOTTOM SUBSTRATE
    n4_power = 0
    y_cur = y_end
    while y_cur < divisions:
        n4_power += middle_power_x * y_field_vals[y_cur]**2
        y_cur += 1
    n4_power *= dy

    total_power = n1_power + n2_power + n3_power + n4_power + n5_power

    x_start_1 = 0
    while x_vals[x_start_1] <= -1 - roughness / w:
        x_start_1 += 1
    x_end_1 = x_start_1
    while x_vals[x_end_1] <= -1 + roughness / w:
        x_end_1 += 1

    y_start_1 = 0
    while y_vals[y_start_1] <= -1 - roughness / h:
        y_start_1 += 1
    y_end_1 = y_start_1
    while y_vals[y_end_1] <= -1 + roughness / h:
        y_end_1 += 1

    y_clad = 0
    cladding_1 = 0
    y_cur = y_start_1
    while y_cur < y_end_1:
        y_clad += y_field_vals[y_cur]**2
        y_cur += 1
    y_clad *= dy

    x_cur = x_start
    while x_cur < x_end:
        cladding_1 += x_field_vals[x_cur]**2 * y_clad
        x_cur += 1
    cladding_1 *= dx

    x_clad = 0
    x_cur = x_start_1
    while x_cur < x_end_1:
        x_clad += x_field_vals[x_cur]**2
        x_cur += 1
    x_clad *= dx

    y_cur = y_start
    cladding_2 = 0
    while y_cur < y_end:
        y_clad += y_field_vals[y_cur]**2 * x_clad
        y_cur += 1
    cladding_2 *= dy

    return 2 * (cladding_1 + cladding_2) / total_power

def confinement_factor_2d(lamb=wavelength, mode=0, mode_2=0, print_statements=False):
    """Computes the confinement factor looking a single mode

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    divisions : number
        number of points to evaluate
    print_statements : bool
        whether or not to print at certain steps of calculation

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a strip waveguide

    """
    res = solve_2d_default(lamb=lamb, mode=mode, mode_2=mode_2, ret_trans=True)
    if res == None:
        return None
    beta, kx, ky = res

    film_x = func_2d(0, kx) - func_2d(-w, kx)
    film_y = func_2d(0, ky) - func_2d(-h, ky)
    total_x = film_x + 1 / ys_lamb_2d_x(kx, wavelength)
    total_y = film_y + 1 / ys_lamb_2d_y(ky, wavelength)
    total = total_x * total_y
    film = film_x * film_y
    if print_statements:
        print("film_power", film, "total_power", total)
    return film / total

def func_2d(x, kx):
    """Helper function for confinement factor calculation

    NOTES
    -----
    • used to perform an integral

    """
    gam = ys_lamb_2d_x(kx, wavelength)
    kappa = kx

    part2 = x / 2 * (gam**2 / kappa**2 + 1)
    part3 = gam / (2 * kappa**2) * cos(2 * kappa * x)
    part4 = sin(2 * kappa * x) / (4 * kappa) * (1 - (gam / kappa)**2)
    return part2 + part3 + part3

def total_confinement_factor(lamb=wavelength, margin=50):
    """Computes the confinement factor looking at all modes

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    margin : number
        x_vals will be within the range     margin_value * w ± w

    RETURNS
    -------
    confinement factor : number
    
    NOTES
    -----
    • Not completely accurate as it does so before performing eigenmode expansion

    """
    freq = c_light / lamb
    omega = freq * 2 * pi
    mode = 0
    res = ey_default(lamb=lamb, mode=mode, margin=margin)
    if res == None:
        return None
    total = 0
    film = 0
    while res != None:
        x_vals, field_vals, beta = res
        x_vals = [i * w for i in x_vals]
        coeff = beta / (2 * omega * u0)
        ec = 1
        phase = phi(beta, lamb)
        ef = ec / cos(phase)
        es = ef * cos(phase - w * kf_lamb(beta, lamb))
        f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + w))
        total += integrate(x_vals, end=-w, func=f1)
        f2 = lambda x: coeff * ef**2 * (cos(kf_lamb(beta, lamb) * x + phase))**2
        temp = integrate(x_vals, start=-w, end=0, func=f2)
        film += temp
        total += temp
        f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
        total += integrate(x_vals, start=0, func=f3)
        mode += 1
        res = ey(mode=mode, margin=20)
    return film / total

"""
------------------------------------------------------------------------------------------------------------
UTILITY FUNCTIONS
------------------------------------------------------------------------------------------------------------
"""

def make_range(start, end, items, log_opt=False):
    """Makes a list of values from start to end (inclusive) with number of items specified

    PARAMETERS
    ----------
    start : number
        first value
    end : number
        last value
    items : int
        number of items in the range
    log_opt : bool
        whether or not to have logarithmic distribution of numbers (as opposed to the default linear)

    RETURNS
    -------
    list of values

    """
    items = (int) (items)
    if log_opt:
        start = log(start, 10)
        end = log(end, 10)
        items = (end - start) / (items - 1)
        toReturn = np.arange(start, end + 0.9*items, items)
        return [10**i for i in toReturn]
    items = (end - start) / (items - 1)
    return list(np.arange(start, end + 0.9*items, items))

def arange(start, end, increment, inclusive=True):
    """Makes a list of values from start to end with increment specified

    PARAMETERS
    ----------
    start : number
        first value
    end : number
        last value
    increment : int
        number of items in the range
    inclusive : bool
        whether or not the list be inclusive of the end boundary

    RETURNS
    -------
    list of values

    """
    if inclusive:
        return list(np.arange(start, end + increment * 0.5, increment))
    else:
        return list(np.arange(start, end, increment))

def lamb_from_V(V_val):
    """Calculates the wavelength corresponding to some wavelength

    PARAMETERS
    ----------
    V_val : number
        V

    RETURNS
    -------
    lambda : number
        wavelength

    """
    return (2 * pi * w * sqrt(nf**2 - ns**2)) / V_val

wavelength_from_V = lamb_from_V
w_from_V = lamb_from_V

def find_intersection(lst1, lst2):
    """Finds the index at which the lists intersect

    PARAMETERS
    ----------
    lst1 : list of numbers

    lst2 : list of numbers

    RETURNS
    -------
    index : int
        index at which they intersect

    NOTES
    -----
    • The lists can only intersect once
    • If the lists do not intersect, the function returns None

    """
    t = 0
    last = len(lst1) - 1
    if (lst1[0] > lst2[0] and lst1[last] > lst2[last]) or (lst1[0] < lst2[0] and lst1[last] < lst2[last]) or lst1[0] == lst1[last] or lst2[0] == lst2[last]:
        return # return None
    while lst1[t] != lst2[t]:
        if lst1[t] > lst2[t]:
            # lst1 is bigger
            for i in range(len(lst1)):
                if lst1[i] < lst2[i]:
                    return i
        else:
            for i in range(len(lst1)):
                if lst1[i] > lst2[i]:
                    return i

def remove_none(lst1, lst2):
    """looks for all None values present in the SECOND list and removes them from both

    PARAMETERS
    ----------
    lst1 : list

    lst2 : list

    RETURNS
    -------
    toReturn1 : list
        lst1 with values corresponding to None removed
    toReturn2 : list
        lst2 with None values removed

    """
    toReturn1, toReturn2 = [], []
    for i in range(len(lst2)):
        if lst2[i] != None:
            toReturn1 += [lst1[i]]
            toReturn2 += [lst2[i]]
    return toReturn1, toReturn2

def curve_area(x_vals, y_vals, start=-inf, end=inf, total_area=True):
    """computes the area under some curve y = f(x)

    PARAMETERS
    ----------
    x_vals : list of numbers
        list of x values
    y_vals : list of numbers
        list of y values
    start : number
        where to start evaluating area
    end : number
        where to stop evaluating area
    total_area : bool
        whether or not to calculate the total area under the curve (take absolute value of points)

    RETURNS
    -------
    area : number
        evaluated area contained by the curve

    """
    x_start = x_vals[0]
    x_end = x_vals[len(x_vals) - 1]
    delta = (x_end - x_start) / len(x_vals)
    new_y = []
    if total_area == False:
        func = lambda x:x
    else:
        func = abs
    new_y = 0
    for i in range(len(x_vals)):
        x = x_vals[i]
        y = y_vals[i]
        if x >= start and x <= end:
            new_y += func(y)
    return new_y * delta

def integrate(x_vals, start=-inf, end=inf, func=lambda x:x):
    """computes the approximated integral of a curve
    
    PARAMETERS
    ----------
    x_vals : list of numbers
        list of x values
    start : number
        where to start evaluating area
    end : number
        where to stop evaluating area
    func : function
        some function f(x) to compute the integral of
    
    RETURNS
    -------
    integral : number
        approximated integral

    """
    x_start = x_vals[0]
    x_end = x_vals[len(x_vals) - 1]
    delta = (x_end - x_start) / len(x_vals)
    output = 0
    for i in range(len(x_vals)):
        x = x_vals[i]
        if x >= start and x <= end:
            output += func(x)
    return output * delta

def optical_power(lamb=wavelength, mode=0, opt=5):
    """Returns the power in one section of the film
    
    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    opt : int
        option
        OPT     REGION
        1       substrate
        2       film
        3       cover
        4       substrate + cover
        5       substrate + film + cover

    RETURNS
    -------
    power : number
        optical power for the given option

    NOTES
    -----
    • If no beta value can be evaluated, it returns None
    • Assumes a slab waveguide

    """
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_default(lamb=wavelength, mode=mode, margin=20)
    if res == None:
        return # return None
    x_vals, field_vals, beta = res
    x_vals = [i * w for i in x_vals]
    coeff = beta / (2 * omega * u0)
    ec = 1
    phase = phi(beta, lamb)
    ef = ec / cos(phase)
    es = ef * cos(phase - w * kf_lamb(beta, lamb))
    total = 0
    if opt == 1:
        f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + w))
        return integrate(x_vals, end=-w, func=f1)
    if opt == 2:
        f2 = lambda x: coeff * ef**2 * (cos(kf_lamb(beta, lamb) * x + phase))**2
        return integrate(x_vals, start=-w, end=0, func=f2)
    if opt == 3:
        f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
        return integrate(x_vals, start=0, func=f3)
    if opt == 4:
        f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + w))
        f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
        total = integrate(x_vals, end=-w, func=f1)
        return total + integrate(x_vals, start=0, func=f3)
    if opt == 5:
        f1 = lambda x: coeff * es**2 * exp(2 * ys_lamb(beta, lamb) * (x + w))
        f2 = lambda x: coeff * ef**2 * (cos(kf_lamb(beta, lamb) * x + phase))**2
        f3 = lambda x: coeff * ec**2 * e**(-2 * ys_lamb(beta, lamb)* x)
        total = integrate(x_vals, end=-w, func=f1) 
        total += integrate(x_vals, start=-w, end=0, func=f2)
        return total + integrate(x_vals, start=0, func=f3)

def cutoff_V(mode=0):
    """Calculates the cutoff V value for some mode

    PARAMETERS
    ----------
    mode : int
        mode number

    RETURNS
    -------
    V : number
        cutoff V value

    """
    return mode * pi + atan((ns**2 - nc**2)/(nf**2 / ns**2))

def cutoff_nf(mode=0):
    """Calculates the cutoff film refractive index value for some mode

    PARAMETERS
    ----------
    mode : int
        mode number

    RETURNS
    -------
    nf : number
        cutoff nf (film refractive index) value

    """
    return sqrt((cutoff_V(mode=mode) / (k * w))**2 + ns**2)

def save_to_file(target_file):
    """Saves the plotted image to a file for display

    PARAMETERS
    ----------
    target_file : string
        target file name to save the image

    """
    savefig(target_file)

"""
------------------------------------------------------------------------------------------------------------
"""