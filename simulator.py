# IMPORT STATEMENTS...take some out
from waveguide import *
import waveguide
from utils import *
from parameters import *
from solvers import *
from fields import *
from bV import *
from loss import *

"""
------------------------------------------------------------------------------------------------------------
SETTING WAVEGUIDE PARAMETERS
------------------------------------------------------------------------------------------------------------
"""

def set_n_film_mat(new_film):
    """Set the film refractive index

    PARAMETERS
    ----------
    new_nf : number
        desired film refractive index

    """
    waveguide.n_film_mat = new_film
    waveguide.nf = refrac[waveguide.method - 1](waveguide.duty_cycle)

def set_n_sub_mat(new_ns):
    """Set the substrate refractive index

    PARAMETERS
    ----------
    new_ns : number
        desired substrate refractive index

    """
    waveguide.n_sub_mat = new_ns
    waveguide.ns=waveguide.nc=waveguide.n2=waveguide.n3=waveguide.n4=waveguide.n5 = new_ns

def set_method(method):
    """
    """
    if method > 3 or method < 1:
        return
    waveguide.method = method
    waveguide.nf=waveguide.n1=refrac[method - 1](waveguide.duty_cycle)

def set_duty_cycle(new_duty_cycle):
    """
    """
    if new_duty_cycle > 1 or new_duty_cycle < 0:
        return
    waveguide.duty_cycle = new_duty_cycle
    waveguide.nf = waveguide.n1 = refrac[waveguide.method - 1](new_duty_cycle)


def set_width(new_w):
    """Set the width of the waveguide to simulate

    PARAMETERS
    ----------
    new_w : number
        desired waveguide width

    """
    waveguide.w=waveguide.default_w=waveguide.fixed_w = new_w

def set_height(new_h):
    """Set the width of the waveguide to simulate

    PARAMETERS
    ----------
    new_h : number
        desired waveguide height

    """
    waveguide.h=waveguide.default_h=waveguide.fixed_h = new_h

def set_doping(doping):
    if doping != "n" and doping != "p":
        return
    else:
        waveguide.doping = doping

def set_concentration(N):
    if N < 0:
        return
    else:
        waveguide.concentration = N

def set_roughness(roughness):
    waveguide.roughness = roughness

def set_wavelength(wavelength):
    waveguide.waveguide = wavelength

def restore_values():
    """Restore waveguide parameters to their default values

    NOTES
    -----
    • Default width : 450nm
    • Default height : 200nm
    • Default nf : 2.48
        - equivalent to refrac2(0.4) with polysi film and sio2 substrate, an swg waveguide with duty cycle of 0.4

    """
    waveguide.default_w=waveguide.w=waveguide.fixed_w = 450e-9
    waveguide.default_h=waveguide.h=waveguide.fixed_h = 200e-9
    waveguide.ns=waveguide.nc=waveguide.n2=waveguide.n3=waveguide.n4=waveguide.n5 = sio2
    waveguide.n_film_mat = si
    waveguide.nf=waveguide.n1 = refrac2(0.4)
    waveguide.method = 2
    waveguide.duty_cycle = 0.4
    waveguide.doping = "n"
    waveguide.concentration = 1e18
    display_variables()

def display_variables():
    """Prints the values of the current waveguide parameters

    """
    print("w:", waveguide.w)
    print("h:", waveguide.h)
    print("nf:", waveguide.nf)
    print("ns:", waveguide.ns)
    print("doping:", waveguide.doping, "at", waveguide.concentration)
    print("duty cycle:", waveguide.duty_cycle)

"""
------------------------------------------------------------------------------------------------------------
"""