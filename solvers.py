from waveguide import *
import waveguide

# CONSTANT SOLVER PARAMETERS
default_iter=10
divisions=500
plot_divisions=5000
print_N=False
print_neff=False

"""
------------------------------------------------------------------------------------------------------------
SOLVERS
------------------------------------------------------------------------------------------------------------
"""
# SOLVER DIMENSION DESCRIBES THE NUMBER OF NON-UNIFORM AXES OF THE WAVEGUIDE CROSS-SECTION

# f_oddNumbers: left hand side         f_evenNumbers: right hand side

# (DEFAULT)
def f1(x):
    return waveguide.w * kf(x)

def f1_lamb(x, lamb):
    return waveguide.w * kf_lamb(x, lamb)

def f1_lamb_alt(x, lamb):
    return waveguide.h * kf_lamb(x, lamb)


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


def solve_1d_x(lamb=None, mode=0, print_N=print_N, max_iter=default_iter):
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
    if lamb == None:
        lamb = waveguide.wavelength
    start, end = k_lamb(lamb) * waveguide.ns + 10, k_lamb(lamb) * waveguide.nf - 10
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

def solve_1d_y(lamb=None, mode=0, print_N=print_N, max_iter=default_iter):
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
    if lamb == None:
        lamb = waveguide.wavelength
    start, end = k_lamb(lamb) * waveguide.ns + 10, k_lamb(lamb) * waveguide.nf - 10
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

# (SOLVING FOR TM MODE)
f1_mag = f1_lamb

def f2_mag(x, lamb, mode):
    kappa_f = kf_lamb(x, lamb)
    return atan((waveguide.nf / waveguide.nc)**2 * yc_lamb(x, lamb) / kappa_f) + atan((waveguide.nf / waveguide.ns)**2 * ys_lamb(x, lamb) / kappa_f) + pi * mode

def solve_tm(lamb=None, mode=0, print_N=print_N, max_iter=default_iter):
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
    • Calculates as a slab waveugide along the x axis (uniform along the y axis)

    """
    if lamb == None:
        lamb = waveguide.wavelength
    start, end = k_lamb(lamb) * waveguide.ns + 10, k_lamb(lamb) * waveguide.nf - 10
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
    return waveguide.w * k_lamb(lamb) * sqrt(waveguide.n1**2 - neff**2)

def f1_eff_alt(neff, lamb):
    return waveguide.h * k_lamb(lamb) * sqrt(waveguide.n1**2 - neff**2)

def f3_eff(N, neff, lamb):
    return waveguide.h * k_lamb(lamb) * sqrt(neff**2 - N**2)

def f3_eff_alt(N, neff, lamb):
    return waveguide.w * k_lamb(lamb) * sqrt(neff**2 - N**2)


def f2_eff(neff, lamb, mode):
    return mode * pi + atan(sqrt(neff**2 - waveguide.n3**2) / sqrt(waveguide.n1**2 - neff**2)) + atan(sqrt(neff**2 - waveguide.n5**2) / sqrt(waveguide.n1**2 - neff**2))

def f2_eff_alt(neff, lamb, mode):
    return mode * pi + atan(sqrt(neff**2 - waveguide.n2**2) / sqrt(waveguide.n1**2 - neff**2)) + atan(sqrt(neff**2 - waveguide.n4**2) / sqrt(waveguide.n1**2 - neff**2))

def f4_eff(N, neff, lamb, mode):
    return mode * pi + atan(sqrt(N**2 - waveguide.n2**2) / sqrt(neff**2 - N**2)) + atan(sqrt(N**2 - waveguide.n4**2) / sqrt(neff**2 - N**2))

def f4_eff_alt(N, neff, lamb, mode):
    return mode * pi + atan(sqrt(N**2 - waveguide.n3**2) / sqrt(neff**2 - N**2)) + atan(sqrt(N**2 - waveguide.n5**2) / sqrt(neff**2 - N**2))


def solve_eff(lamb=None, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_neff=print_neff, print_N=print_N, max_iter=default_iter):
    """Solves for Beta using the effective index method (considering a rectangular waveguide)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    if lamb == None:
        lamb = waveguide.wavelength
    start, end = max(waveguide.n2 + 1e-3, waveguide.n4 + 1e-3), waveguide.n1 - 1e-3
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

def solve_eff_alt(lamb=None, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_neff=print_neff, print_N=print_N, max_iter=default_iter):
    """Solves for Beta using the effective index method (considering a rectangular waveguide)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    if lamb == None:
        lamb = waveguide.wavelength
    start, end = max(waveguide.n3 + 1e-3, waveguide.n5 + 1e-3), waveguide.n1 - 1e-3
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

    start, end = max(waveguide.n2 + 1e-3, waveguide.n4 + 1e-3), neff - 1e-3
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
    return waveguide.w * kx

def f7(ky):
    return waveguide.h * ky


def f6(kx, lamb, mode):
    return mode * pi + atan(sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.n2**2) - kx**2) / kx) + atan(sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.n4**2) - kx**2) / kx)

def f8(ky, lamb, mode):
    return mode * pi + atan(sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.n3**2) - ky**2) / ky) + atan(sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.n5**2) - ky**2) / ky)


def solve_marcatili(lamb=None, mode=0, mode_2=0, ret_N=False, ret_trans=False, print_N=print_N, print_trans=False, max_iter=default_iter, ignore_beta=False):
    """Solves for Beta using the Marcatili method (considering a rectangular waveguide)
       
    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    if lamb == None:
        lamb = waveguide.wavelength
    f_right_1 = f6
    f_right_2 = f8
    start, end = 10, sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.ns**2)) - 10
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

    start, end = 10, sqrt(k_lamb(lamb)**2 * (waveguide.n1**2 - waveguide.ns**2)) - 10
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

    res = (k_lamb(lamb) * waveguide.n1)**2 - kx**2 - ky**2
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

def solve_2d(lamb=None, mode=0, mode_2=0, ret_trans=False, max_iter=default_iter, print_statements=False):
    """Modified solver of choice to solve for Beta

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    if lamb == None:
        lamb = waveguide.wavelength
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

    temp = sqrt(((k_lamb(lamb) * waveguide.nf)**2 - (beta)**2) / (kx**2 + ky**2))
    kx *= temp
    ky *= temp

    if ret_trans:
        return beta, kx, ky
    return beta

