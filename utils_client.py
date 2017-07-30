from utils import *
import utils

def loss(duty_cycle, roughness=11e-9, N=1e18, doping="n", constant_absorption=False, lamb=1550e-9, opt=0, factor_fundamental=True, return_components=False, roughness_opt=0, alter_thickness=True, alter_confinement=True, correct_theta=False, index_opt=2, exp=1.813):
    """Calculates the expected loss of a swg waveguide
    
    PARAMETERS
    ----------
    duty cycle : number

    roughness : number
        roughness variance
    N : number
        doping concentration (dopants/cm^2)
    doping : 'n' or 'p'
        doping type
    constant_absorption : bool
        whether or not to have absorption loss vary independent of the film refractive index
    lamb : number
        wavelength (lambda)
    opt : int
        option
        OPT     LOSS MECHANISMS
        0       ROUGHNESS + ABSORPTION
        1       ROUGHNESS
        2       ABSORPTION
    factor_fundamental : bool
        whether or not to perform eigenmode expansion
        (considering all modes besides the fundamental as completely loss)
    return_components : bool
        whether or not to return the individual components of loss
    roughness_opt : opt
        OPT     METHOD
        0       maintain constant roughness
        1       vary roughness directly with duty cycle
        2       vary roughness with the square root of the duty cycle
    alter_thickness : bool
        whether or not to factor SWG waveguide expansion in beta calculation
    alter_confinement : bool
        whether or not to factor SWG waveguide expansion in calculating confinement
    correct_theta : bool
        whether to use a larger propagation angle (treating it similar to a slab waveguide)
    index_opt : opt
        how to estimate SWG waveguide film refractive index (modeling as a strip waveguide)
        OPT     METHOD
        1       refrac1
        2       refrac2
        3       refrac3
    exp : num
        constant describing how the waveguide expands with the applied SWG duty cycle

    RETURNS
    -------
    loss : number

    NOTES
    -----
    • The main function of interest, combining all simulated components to predict total transmission loss
    • If no valid of film refractive index is chosen, it returns None
    • If it is neither 'n' nor 'p' type doped, it returns None

    """
    if doping != "n" and doping != "p":
        return None

    # saving original values
    original_n = utils.n1

    utils.w, utils.h = utils.default_w, utils.default_h
    utils.fixed_w, utils.fixed_h = utils.default_w, utils.default_h
    
    if index_opt == 1:
        n = refrac1(duty_cycle)
    elif index_opt == 2:
        n = refrac2(duty_cycle)
    elif index_opt == 3:
        n = refrac3(duty_cycle)
    else:
        return
    utils.nf = utils.n1 = n
    if alter_thickness:
        calc_w, calc_h = effective_film_thickness_2d(lamb=lamb, mode=0, mode_2=0, fixed_thickness=not alter_thickness)
        add_w, add_h = (1 - duty_cycle)**exp * (calc_w - utils.w), (1 - duty_cycle)**exp * (calc_h - utils.h)
        utils.w, utils.h = utils.w + add_w, utils.h + add_h
        utils.fixed_w, utils.fixed_h = calc_w, calc_h
    else:
        utils.w, utils.h = effective_film_thickness_2d(lamb=lamb, mode=0, mode_2=0, fixed_thickness=not alter_thickness)
    if roughness_opt == 1:
        roughness *= duty_cycle
    elif roughness_opt == 2:
        roughness *= sqrt(duty_cycle)
    if factor_fundamental:
        r_loss_x, r_loss_y = roughness_loss_2d(lamb=lamb, roughness=roughness, method=2, fixed_thickness=alter_thickness, correct_theta=correct_theta)
        a_loss = absorption_loss(N=N, lamb=lamb, doping=doping, constant=constant_absorption)
        r_loss = sqrt((r_loss_x / fundamental_mode_weight(n))**2 + r_loss_y**2)
    else:
        r_loss, a_loss = roughness_loss_2d(lamb=lamb, roughness=roughness, fixed_thickness=alter_thickness, correct_theta=correct_theta), absorption_loss(N=N, lamb=lamb, doping=doping, constant=constant_absorption)
    if alter_confinement and alter_thickness:
        pass
    elif alter_confinement:
        calc_w, calc_h = effective_film_thickness_2d(lamb=lamb, mode=0, fixed_thickness=not alter_thickness)
        add_w, add_h = (1 - duty_cycle)**exp * (calc_w - utils.w), (1 - duty_cycle)**exp * (calc_h - utils.h)
        utils.w, utils.h = utils.w + add_w, utils.h + add_h
    else:
        utils.w, utils.h = utils.default_w, utils.default_h
    c_factor = confinement_factor_2d(lamb=lamb)

    # restore values
    utils.w, utils.h = utils.default_w, utils.default_h
    utils.fixed_w, utils.fixed_h = utils.default_w, utils.default_h
    utils.nf=utils.n1=original_n
    if return_components:
        return duty_cycle, c_factor, r_loss, a_loss
    if roughness_opt != 0:
        return c_factor * (r_loss + duty_cycle * a_loss)
    return duty_cycle * c_factor * (r_loss + a_loss)


def fundamental_mode_weight(n, n_ref=refrac2(0.4), margin=20):
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
    original_n = utils.nf
    utils.nf = utils.n1 = refrac2(duty_cycle)
    res = ey_default(margin=margin, mode=0)
    if res == None:
        return
    temp = res[1]
    temp = [i**2 for i in temp]
    a = curve_area(res[0], temp)
    temp = [i / a for i in temp]
    m1 = max(temp)

    utils.nf=utils.n1 = n
    res = ey_default(margin=margin, mode=0)
    if res == None:
        return
    temp = res[1]
    temp = [i**2 for i in temp]
    a = curve_area(res[0], temp)
    temp = [i / a for i in temp]
    m2 = max(temp)

    utils.nf = utils.n1 = original_n

    return m1 / m2

def adjusted_roughness_loss(duty_cycle, roughness=11e-9, lamb=1550e-9, factor_fundamental=False, return_components=False, roughness_opt=0, correct_theta=False, index_opt=2, exp=1.813):
    """Calculates the loss due to surface roughness assuming a 2D SWG waveguide

    
    PARAMETERS
    ----------
    duty cycle : number

    roughness : number
        roughness variance
    lamb : number
        wavelength (lambda)
    factor_fundamental : bool
        whether or not to perform eigenmode expansion
        (considering all modes besides the fundamental as completely loss)
    return_components : bool
        whether or not to return the individual components of loss
    roughness_opt : opt
        OPT     METHOD
        0       maintain constant roughness
        1       vary roughness directly with duty cycle
        2       vary roughness with the square root of the duty cycle
    correct_theta : bool
        whether to use a larger propagation angle (treating it similar to a slab waveguide)
    index_opt : opt
        how to estimate SWG waveguide film refractive index (modeling as a strip waveguide)
        OPT     METHOD
        1       refrac1
        2       refrac2
        3       refrac3
    exp : num
        constant describing how the waveguide expands with the applied SWG duty cycle

    Returns
    -------
    loss : number
        loss due to surface roughness and consequent scattering

    NOTES
    -----
    • If no Beta can be evaluated for this mode, it returns None
    • Assumes a strip SWG waveguide with corresponding widened dimensions

    """
    utils.w, utils.h = utils.default_w, utils.default_h
    utils.fixed_w, utils.fixed_h = utils.default_w, utils.default_h
    if index_opt == 1:
        n = refrac1(duty_cycle)
    elif index_opt == 2:
        n = refrac2(duty_cycle)
    elif index_opt == 3:
        n = refrac3(duty_cycle)
    else:
        return
    utils.nf = utils.n1 = n
    calc_w, calc_h = effective_film_thickness_2d(lamb=lamb, mode=0, mode_2=0, fixed_thickness=False)
    add_w, add_h = (1 - duty_cycle)**exp * (calc_w - utils.w), (1 - duty_cycle)**exp * (calc_h - utils.h)
    utils.w, utils.h = utils.w + add_w, utils.h + add_h
    utils.fixed_w, utils.fixed_h = calc_w, calc_h
    if roughness_opt == 1:
        roughness *= duty_cycle
    elif roughness_opt == 2:
        roughness *= sqrt(duty_cycle)
    if factor_fundamental:
        r_loss_x, r_loss_y = roughness_loss_2d(lamb=lamb, roughness=roughness, method=2, fixed_thickness=True, correct_theta=correct_theta)
        r_loss = sqrt((r_loss_x / fundamental_mode_weight(n))**2 + r_loss_y**2)
    else:
        r_loss = roughness_loss_2d(lamb=lamb, roughness=roughness, fixed_thickness=True, correct_theta=correct_theta)
    utils.w, utils.h = utils.default_w, utils.default_h
    utils.fixed_w, utils.fixed_h = utils.default_w, utils.default_h
    return r_loss

"""
------------------------------------------------------------------------------------------------------------
SETTING WAVEGUIDE PARAMETERS
------------------------------------------------------------------------------------------------------------
"""

def set_nf(new_nf):
    """Set the film refractive index

    PARAMETERS
    ----------
    new_nf : number
        desired film refractive index

    """
    utils.nf=utils.n1 = new_nf

def set_ns(new_ns):
    """Set the substrate refractive index

    PARAMETERS
    ----------
    new_ns : number
        desired substrate refractive index

    """
    utils.ns=utils.nc=utils.n2=utils.n3=utils.n4=utils.n5 = new_ns

def set_width(new_w):
    """Set the width of the waveguide to simulate

    PARAMETERS
    ----------
    new_w : number
        desired waveguide width

    """
    utils.w=utils.default_w=utils.fixed_w = new_w

def set_height(new_h):
    """Set the width of the waveguide to simulate

    PARAMETERS
    ----------
    new_h : number
        desired waveguide height

    """
    utils.h=utils.default_h=utils.fixed_h = new_h

def restore_values():
    """Restore waveguide parameters to their default values

    NOTES
    -----
    • Default width : 450nm
    • Default height : 200nm
    • Default nf : 2.48     equivalent to refrac2(0.4), an swg waveguide with duty cycle of 0.4
    """
    utils.default_w=utils.w=utils.fixed_w = 450e-9
    utils.default_h=utils.h=utils.fixed_h = 200e-9
    utils.nf=utils.n1 = refrac2(0.4)
    utils.ns=utils.nc=utils.n2=utils.n3=utils.n4=utils.n5 = utils.sio2
    display_variables()

def display_variables():
    """Prints the values of the current waveguide parameters

    """
    print("w:", utils.default_w)
    print("h:", utils.default_h)
    print("nf:", utils.nf)
    print("ns:", utils.ns)

"""
------------------------------------------------------------------------------------------------------------
"""