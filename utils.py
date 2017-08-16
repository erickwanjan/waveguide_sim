from numpy import *
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import *

"""
------------------------------------------------------------------------------------------------------------
UTILITY FUNCTIONS
------------------------------------------------------------------------------------------------------------
"""

# MATH FUNCTIONS
rad = lambda x: x * pi / 180
scalar = lambda x: x * 180 / pi

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

def save_to_file(target_file):
    """Saves the plotted image to a file for display

    PARAMETERS
    ----------
    target_file : string
        target file name to save the image

    """
    savefig(target_file)
