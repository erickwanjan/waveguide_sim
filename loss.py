import waveguide
from waveguide import *
import parameters
from parameters import *
import solvers
from solvers import *


"""
------------------------------------------------------------------------------------------------------------
LOSS FUNCTIONS
------------------------------------------------------------------------------------------------------------
"""

def absorption_loss(print_afc=False, constant_absorption=False):
    """Calculates the loss due to free carrier absorption
    PARAMETERS
    ----------
    print_afc : bool
        whether or not to print per unit dopant concentration
    constant_absorption : bool
        whether or not loss is kept invariant with the effective index

    RETURNS
    -------
    loss : number
        Absorption loss due to free carrier (intraband) absorption

    NOTES
    -----
    • Assumes that the material waveguide used is heavily doped silicon

    """
    lamb = waveguide.wavelength
    N, doping = waveguide.concentration, waveguide.doping
    
    if constant_absorption == True:
        nf = wavelength.si
    else:
        nf = waveguide.nf
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
        print("afc", afc)
    return afc * N

def roughness_loss_uniform_slab(lamb=None, mode=0, print_statements=False):
    """Calculates the loss due to surface roughness

    Parameters
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    print_statements : bool
        whether or not to print at certain steps of loss calculation

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a uniform slab waveguide along the x-axis (uniform along the y-axis)

    """
    if lamb == None:
        lamb = waveguide.wavelength
    res = propagation_angle(lamb=lamb, mode=mode, ret_beta=True)
    if res == None:
        return
    theta_incident, beta = res
    if theta_incident == None:
        return
    theta_m = 90 - theta_incident
    A = 4 * pi * waveguide.roughness / lamb
    if print_statements:
        print("theta", theta_incident)
        print("A", A)

    return A**2 * (cos(rad(theta_m)))**3 / (2 * sin(rad(theta_m)) * effective_film_thickness_1d_x(beta, lamb)) * 4.343 / 100

def roughness_loss_1d_x(lamb=None, mode=0, return_components=False):
    """Calculates the loss due to surface roughness assuming a slab SWG waveguide uniform in the x direction

    
    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    return_components : bool
        whether or not to return the individual components of loss

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a slab SWG waveguide with corresponding widened dimensions

    """
    roughness = waveguide.roughness
    duty_cycle = waveguide.duty_cycle
    exp = waveguide.adjustment_exp

    if lamb == None:
        lamb = waveguide.wavelength

    utils.w = utils.default_w
    utils.fixed_w = utils.default_w

    calc_w = effective_film_thickness_1d(lamb=lamb, mode=mode)
    add_w = (1 - duty_cycle)**exp * (calc_w - waveguide.w)
    waveguide.w = waveguide.w + add_w
    waveguide.fixed_w, waveguide.fixed_h = calc_w, calc_h

    r_loss = roughness_loss_uniform_slab(lamb=lamb, mode=mode)

    # restore values
    waveguide.w = waveguide.default_w
    waveguide.fixed_w = waveguide.default_w

    return r_loss

def roughness_loss_uniform_strip(lamb=None, mode=0, mode_2=0, fixed_thickness=False, print_statements=False, return_components=False):
    # CHANGE METHOD TO RETURN COMPONENTS
    """Calculates the loss due to surface roughness assuming a 2D waveguide for a specific propagation mode

    Parameters
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    fixed_thickness : bool
        whether or not to assume the penetration depth is 0
    print_statements : bool
        whether or not to print at certain steps of loss calculation
    return_components : bool
        whether or not to return the x and y components of loss

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering
    loss_x : number
        loss due to to surface roughness with the boundary along the x axis
    loss_y : number
        loss due to to surface roughness with the boundary along the y axis

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip waveguide

    """
    roughness = waveguide.roughness
    nf = waveguide.nf

    if lamb == None:
        lamb = waveguide.wavelength
    res = solve_marcatili(lamb=lamb, mode=mode, mode_2=mode_2, ret_trans=True)
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

    theta_x = scalar(atan(kx/beta))
    theta_y = scalar(atan(ky/beta))
    theta_m_x = 90 - theta_x
    theta_m_y = 90 - theta_y
    if fixed_thickness:
        thickness_x = waveguide.fixed_w
        thickness_y = waveguide.fixed_h
    else:
        thickness_x = effective_film_thickness_2d_x(kx, lamb)
        thickness_y = effective_film_thickness_2d_y(ky, lamb)

    loss_x = (cos(rad(theta_m_x)))**3 / (2 * sin(rad(theta_m_x)) * thickness_x)
    loss_y = (cos(rad(theta_m_y)))**3 / (2 * sin(rad(theta_m_y)) * thickness_y)

    if print_statements:
        print("beta:", beta, "kx:", kx, "ky", ky)
        print("theta_x", theta_x, "theta_y", theta_y)
        print("theta_m_x:", theta_m_x, "theta_m_y:", theta_m_y)
        print("thickness_x:", thickness_x, "thickness_y:", thickness_y)
        print("loss_components", A**2 * loss_x * 4.343 / 100, A**2 * loss_y * 4.343 / 100)
        print("per unit roughness", A**2 * (loss_x + loss_y) * 4.343 / 100 / roughness**2)
    if not return_components:
        return A**2 * sqrt(loss_x**2 + loss_y**2) * 4.343 / 100
    else:
        return A**2 * loss_x * 4.343 / 100, A**2 * loss_y * 4.343 / 100

def roughness_loss_2d(lamb=None):
    """Calculates the loss due to surface roughness assuming a 2D SWG waveguide

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip SWG waveguide with corresponding widened dimensions

    """
    roughness = waveguide.roughness
    duty_cycle = waveguide.duty_cycle
    exp = waveguide.adjustment_exp

    if lamb == None:
        lamb = waveguide.wavelength

    utils.w, utils.h = utils.default_w, utils.default_h
    utils.fixed_w, utils.fixed_h = utils.default_w, utils.default_h

    calc_w, calc_h = effective_film_thickness_2d(lamb=lamb, mode=0, mode_2=0)
    add_w, add_h = (1 - duty_cycle)**exp * (calc_w - waveguide.w), (1 - duty_cycle)**exp * (calc_h - waveguide.h)
    waveguide.w, waveguide.h = waveguide.w + add_w, waveguide.h + add_h
    waveguide.fixed_w, waveguide.fixed_h = calc_w, calc_h

    r_loss_x, r_loss_y = roughness_loss_uniform_strip(lamb=lamb, mode=0, mode_2=0, fixed_thickness=True, return_components=True)
    r_loss = sqrt((r_loss_x / fundamental_mode_weight(n))**2 + r_loss_y**2)

    # restore values
    waveguide.w, waveguide.h = waveguide.default_w, waveguide.default_h
    waveguide.fixed_w, waveguide.fixed_h = waveguide.default_w, waveguide.default_h

    return r_loss

def transmission_loss(lamb=1550e-9, return_components=False, constant_absorption=False):
    """Calculates the expected loss of a swg waveguide
    
    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for propagation in x direction
    mode_2 : int
        mode number for propagation in y direction
    return_components : bool
        whether or not to return the individual components of loss
    constant_absorption : bool

    RETURNS
    -------
    loss : number

    NOTES
    -----
    • The main function of interest, combining all simulated components to predict total transmission loss
    • If no valid of film refractive index is chosen, it returns None
    • If it is neither 'n' nor 'p' type doped, it returns None

    """
    N, doping = waveguide.N, waveguide.doping
    roughness = waveguide.roughness
    duty_cycle = waveguide.duty_cycle
    exp = waveguide.adjustment_exp

    if lamb == None:
        lamb = waveguide.wavelength

    res = effective_film_thickness_2d(lamb=lamb, mode=0, mode_2=0)
    if res == None:
        return
    calc_w, calc_h = res
    add_w, add_h = (1 - duty_cycle)**exp * (calc_w - waveguide.w), (1 - duty_cycle)**exp * (calc_h - waveguide.h)
    waveguide.w, waveguide.h = waveguide.w + add_w, waveguide.h + add_h
    waveguide.fixed_w, waveguide.fixed_h = calc_w, calc_h

    r_loss_x, r_loss_y = roughness_loss_uniform_strip(lamb=lamb, mode=0, mode_2=0, fixed_thickness=True, return_components=True)
    a_loss = absorption_loss(constant_absorption=constant_absorption)
    c_factor = confinement_factor_2d(lamb=lamb)
    r_loss = sqrt((r_loss_x / fundamental_mode_weight(lamb=lamb))**2 + r_loss_y**2)
    print(waveguide.fixed_w, waveguide.fixed_h, waveguide.w, waveguide.h)
    # restore values
    waveguide.w, waveguide.h = waveguide.default_w, waveguide.default_h
    waveguide.fixed_w, waveguide.fixed_h = waveguide.default_w, waveguide.default_h

    if return_components:
        return duty_cycle * c_factor * r_loss, duty_cycle * c_factor * a_loss
    return duty_cycle * c_factor * (r_loss + a_loss)
