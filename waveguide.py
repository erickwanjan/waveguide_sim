"""
TODO
    Add sample code to run functions (savefig)
    
    waveguide class...
    file for solvers
    file for fields/confinement
    file for angles, film thicknesses
    file for loss
    utils file
    lamb=wavelength thing...
"""
from utils import *

"""
------------------------------------------------------------------------------------------------------------
INITIAL PARAMETERS
------------------------------------------------------------------------------------------------------------
"""
si = 3.5
sio2 = 1.45
air = 1

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

refrac = [refrac1, refrac2, refrac3]

def recover_duty_cycle_1(n):
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
    • If it is not possible to create a waveguide with this refractive index for the given materials, it returns None
    • Works for the approximation used in refrac1

    """
    if n < ns or n > nf:
        return
    return (n - ns) / (nf - ns)

def recover_duty_cycle_2(n):
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
    • If it is not possible to create a waveguide with this refractive index for the given materials, it returns None
    • Works for the approximation used in refrac2

    """
    if n < ns or n > nf:
        return
    return (n**2 - ns**2) / (nf**2 - ns**2)

def recover_duty_cycle_3(n):
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
    • If it is not possible to create a waveguide with this refractive index for the given materials, it returns None
    • Works for the approximation used in refrac3

    """
    if n < ns or n > nf:
        return
    return (sqrt(n) - sqrt(ns)) / (sqrt(nf) - sqrt(ns))

def recover_duty_cycle(n):
    """Estimates the duty cycle of an SWG with the given effective index according to your chosen method

    PARAMETERS
    ----------
    n : number
        effective refractive index

    RETURNS
    -------
    duty cycle : number

    NOTES
    -----
    • If it is not possible to create a waveguide with this refractive index for the given materials, it returns None

    """
    if method == 1:
        return recover_duty_cycle_1(n)
    elif method == 2:
        return recover_duty_cycle_2(n)
    elif method == 3:
        return recover_duty_cycle_3(n)

duty_cycle = 0.4
method = 2

n_film_mat = si
n_sub_mat = sio2

nf = refrac[method - 1](duty_cycle)
nc = n_sub_mat
ns = n_sub_mat

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

adjustment_exp = 1.813

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

concentration = 1e18
doping = "n"
roughness=11e-9

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

kf = lambda beta: sqrt(k**2 * nf**2 - beta**2)
ys = lambda beta: sqrt(beta**2 - k**2 * ns**2)
yc = lambda beta: sqrt(beta**2 - k**2 * nc**2)
kf_lamb = lambda beta, lamb: sqrt(k_lamb(lamb)**2 * nf**2 - beta**2)
ys_lamb = lambda beta, lamb: sqrt(beta**2 - k_lamb(lamb)**2 * ns**2)
yc_lamb = lambda beta, lamb: sqrt(beta**2 - k_lamb(lamb)**2 * nc**2)
ys_lamb_2d_x = lambda kx, lamb: sqrt(k_lamb(lamb)**2 * (n1**2 - n2**2) - kx**2)
ys_lamb_2d_y = lambda ky, lamb: sqrt(k_lamb(lamb)**2 * (n1**2 - n2**2) - ky**2)

phi = lambda beta, lamb: atan(yc_lamb(beta, lamb) / kf_lamb(beta, lamb)) # phase constant

def lamb_from_V(V_val):
    """Calculates the wavelength corresponding to some V Value (in bV curve plotting)

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

"""
------------------------------------------------------------------------------------------------------------
"""