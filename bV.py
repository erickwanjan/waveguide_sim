import solvers
from solvers import *

"""
------------------------------------------------------------------------------------------------------------
BV CURVES
------------------------------------------------------------------------------------------------------------
"""

def bV(start=1e-8, end=1e-5, mode=0, iterations=200, solver=solve_1d_x, max_solver_iter=default_iter):
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
    V = [k_lamb(l) * waveguide.w * sqrt(waveguide.nf**2 - waveguide.ns**2) for l in lamb_lst]
    N_lst = [N(beta_lst[i], lamb_lst[i]) for i in range(len(lamb_lst))]
    b = [(N_lst[i]**2 - waveguide.ns**2) / (waveguide.nf**2 - waveguide.ns**2) for i in range(len(lamb_lst))]
    return V, b

def bV_2d_lamb(start=1e-8, end=1e-5, mode=0, mode_2=0, iterations=200, solver=solve_2d, max_solver_iter=default_iter):
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
    V = [k_lamb(l) * waveguide.w * sqrt(waveguide.nf**2 - waveguide.ns**2) for l in lamb_lst]
    N_lst = [N(beta_lst[i], lamb_lst[i]) for i in range(len(lamb_lst))]
    b = [(N_lst[i]**2 - waveguide.ns**2) / (waveguide.nf**2 - waveguide.ns**2) for i in range(len(lamb_lst))]
    return V, b
