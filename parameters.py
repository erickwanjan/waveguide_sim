from fields import *
import fields

"""
------------------------------------------------------------------------------------------------------------
ANGLES
------------------------------------------------------------------------------------------------------------
"""

def propagation_angle(lamb=None, mode=0, ret_beta=False):
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
    if lamb == None:
        lamb = waveguide.wavelength
    beta = solve_1d_x(lamb, mode)
    if beta == None:
        return None
    theta = scalar(atan(kf_lamb(beta, lamb) / beta))
    if ret_beta:
        return theta, beta
    return theta

def propagation_angle_2d(lamb=None, mode=0, mode_2=0, ret_beta=False):
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
    if lamb == None:
        lamb = waveguide.wavelength
    res = solve_2d(lamb=lamb, mode=mode, ret_trans=True)
    if res == None:
        return None
    beta, kx, ky = res
    if ret_beta:
        return beta, scalar(atan(kx / beta)), scalar(atan(ky / beta))
    return scalar(atan(kx / beta)), scalar(atan(ky / beta))

def critical_angle():
    """Calculates the critical angle of a slab waveguide
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

    NOTES
    -----
    • The critical angle is defined between light and the vector normal to the waveguide boundary

    """
    return scalar(asin(ns / nf))

def critical_angle_incident():
    """Calculates the critical angle of a slab waveguide
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
    
    NOTES
    -----
    • The critical angle is defined between light and the vector normal to the waveguide boundary

    """
    return 90 - critical_angle()

"""
------------------------------------------------------------------------------------------------------------
EFFECTIVE FILM THICKNESS
------------------------------------------------------------------------------------------------------------
"""

def effective_film_thickness_1d_x(lamb=None, solver=solve_1d_x):
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
    if lamb == None:
        lamb = waveguide.wavelength
    beta = solver(mode=mode, lamb=lamb)
    if beta == None:
        return
    return w + 2 / ys_lamb(beta, lamb)

def effective_film_thickness_2d_x(kx, lamb=None):
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
    if lamb == None:
        lamb = waveguide.wavelength
    w, n1, n2 = waveguide.w, waveguide.n1, waveguide.n2
    if kx**2 >= k_lamb(lamb)**2 * (n1**2 - n2**2):
        return w
    return w + 2 / ys_lamb_2d_x(kx, lamb)

def effective_film_thickness_2d_y(ky, lamb=None):
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
    if lamb == None:
        lamb = waveguide.wavelength
    h, n1, n2 = waveguide.h, waveguide.n1, waveguide.n2
    if ky**2 >= k_lamb(lamb)**2 * (n1**2 - n2**2):
        return h
    return h + 2 / ys_lamb_2d_y(ky, lamb)

def effective_film_thickness_2d(lamb=None, mode=0, mode_2=0):
    """computes the effective thickness of the waveguide in the x and y direction

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction

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
    if lamb == None:
        lamb = waveguide.wavelength
    nf = waveguide.nf
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

    thickness_x = effective_film_thickness_2d_x(kx, lamb)
    thickness_y = effective_film_thickness_2d_y(ky, lamb)
    return thickness_x, thickness_y # thickness_x (w), thickness_y (h)

"""
------------------------------------------------------------------------------------------------------------
CONFINEMENT FACTOR
------------------------------------------------------------------------------------------------------------
"""

def confinement_factor_single_mode(lamb=None, solver=solve_1d_x):
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
    if lamb == None:
        lamb = waveguide.wavelength
    beta = solver(lamb=lamb)
    yc_val = yc_lamb(beta, lamb)
    ys_val = ys_lamb(beta, lamb)
    kf_val = kf_lamb(beta, lamb)
    top = w + yc_val / (kf_val**2 + yc_val**2) + ys_val / (kf_val**2 + ys_val**2)
    bot = w + 1 / yc_val + 1 / ys_val
    return top / bot

def confinement_factor_x(lamb=None, mode=0, margin=20, print_statements=False):
    """Computes the confinement factor looking at a single mode of choice

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
    if lamb == None:
        lamb = waveguide.wavelength
    w = waveguide.w

    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_1d_x(lamb=lamb, mode=mode, margin=margin)
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

def confinement_factor_y(lamb=None, mode=0, margin=20, print_statements=False):
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
    if lamb == None:
        lamb = waveguide.wavelength
    h = waveguide.h

    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_1d_y(lamb=lamb, mode=mode, margin=margin)
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

def confinement_factor_interface(roughness=11.5e-9, lamb=None, mode=0, margin=20, print_statements=False):
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
    if lamb == None:
        lamb = waveguide.wavelength
    w = waveguide.w
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_1d_x(lamb=lamb, mode=mode, margin=margin)
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

def confinement_factor_interface_2d(roughness=11.5e-9, lamb=None, mode=0, mode_2=0, margin=20, divisions=None, correct_2d=True):
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

    RETURNS
    -------
    confinement factor : number

    NOTES
    -----
    • If no beta can be evaluated, it returns None
    • Assumes a strip waveguide

    """
    if lamb == None:
        lamb = waveguide.wavelength
    w, h = waveguide.w, waveguide.h
    if divisions == None:
        divisions = solvers.plot_divisions

    divisions = (divisions) * 1 * margin
    freq = c_light / lamb
    omega = freq * 2 * pi
    if mode_2 != 0:
        mod_correct = False
    else:
        mod_correct=correct_2d
    res = ey_1d_x(lamb=lamb, mode=mode, margin=margin, divisions=divisions, correct_2d=mod_correct)
    if res == None:
        return
    x_vals, x_field_vals, beta = res
    if mode_2 != 0:
        mod_correct = False
    else:
        mod_correct=correct_2d
    res = ey_1d_y(lamb=lamb, mode=mode_2, margin=margin, divisions=divisions, correct_2d=mod_correct)
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

def confinement_factor_2d(lamb=None, mode=0, mode_2=0, print_statements=False):
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
    if lamb == None:
        lamb = waveguide.wavelength
    w, h = waveguide.w, waveguide.h
    res = solve_2d(lamb=lamb, mode=mode, mode_2=mode_2, ret_trans=True)
    if res == None:
        return None
    beta, kx, ky = res

    film_x = func_2d(0, kx) - func_2d(-w, kx)
    film_y = func_2d(0, ky) - func_2d(-h, ky)
    total_x = film_x + 1 / ys_lamb_2d_x(kx, lamb)
    total_y = film_y + 1 / ys_lamb_2d_y(ky, lamb)
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

"""
------------------------------------------------------------------------------------------------------------
OPTICAL POWER
------------------------------------------------------------------------------------------------------------
"""

def optical_power(lamb=None, mode=0, opt=5):
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
    • Assumes a slab waveguide along the x axis

    """
    if lamb == None:
        lamb = waveguide.wavelength
    w = waveguide.w
    freq = c_light / lamb
    omega = freq * 2 * pi
    res = ey_1d_x(lamb=wavelength, mode=mode, margin=20)
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

"""
------------------------------------------------------------------------------------------------------------
EIGENMODE EXPANSION
------------------------------------------------------------------------------------------------------------
"""

def fundamental_mode_weight(lamb=None, n_ref=2.4822368944160025, margin=20):
    """Performs eigenmode expansion to determine the fundamental mode weight

    PARAMETERS
    ----------
    n : number
        refractive index of waveguide to perform expansion on
    n_ref : number
        reference refractive index (with only single mode behavior) with fundamental mode fit to Gaussian input wave
    margin : number
        how far from the waveguide film to consider in evaluation

    RETURNS
    -------
    fundamental mode weight : number

    NOTES
    -----
    • If beta cannot be evaluated, it returns None

    """
    if lamb == None:
        lamb = waveguide.wavelength

    original_n = waveguide.nf

    res = ey_1d_x(lamb=lamb, margin=margin, mode=0)
    # res = ey_1d_x(lamb=lamb, margin=margin, mode=0, correct_2d=True)
    if res == None:
        return
    temp = res[1]
    temp = [i**2 for i in temp]
    a = curve_area(res[0], temp)
    temp = [i / a for i in temp]
    m2 = max(temp)

    waveguide.nf=waveguide.n1 = n_ref
    res = ey_1d_x(lamb=lamb, margin=margin, mode=0)
    # res = ey_1d_x(lamb=lamb, margin=margin, mode=0, correct_2d=True)
    if res == None:
        return
    temp = res[1]
    temp = [i**2 for i in temp]
    a = curve_area(res[0], temp)
    temp = [i / a for i in temp]
    m1 = max(temp)

    waveguide.nf=waveguide.n1 = original_n

    return m1 / m2

"""
------------------------------------------------------------------------------------------------------------
"""