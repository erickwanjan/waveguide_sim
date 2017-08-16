from solvers import *
import solvers

def ey_1d_x(lamb=None, mode=0, margin=1.5, max_solver_iter=None, divisions=None, correct_2d=False):
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
        whether or not to produce a more accurate 2d field distribution (though some boundary conditions may not be met)

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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc

    if lamb == None:
        lamb = waveguide.wavelength

    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * w * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_1d_x(lamb=lamb, mode=mode, max_iter=max_solver_iter)
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

def ey_1d_y(lamb=None, mode=0, margin=1.5, max_solver_iter=None, divisions=None, correct_2d=False):
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
        produce a more accurate 2d field distribution (though some boundary conditions may not be met)

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
    • List of y_vals returned are in terms of the height

    """
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * h * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_1d_y(lamb=lamb, mode=mode, max_iter=max_solver_iter)
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

def plot_ey_1d_x(lamb=None, mode=0, margin=1.5, peak_normalize=False, true_normalize=False, center=False, max_solver_iter=None, divisions=None, correct_2d=False, target_file=None, no_output=False):
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
        whether or not to center the plot(moving margins of film to -0.5w, 0.5w)
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    correct_2d : bool
        whether or not to produce a more accurate 2d field distribution (though some boundary conditions may not be met)
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot
    no_output : bool
        whether or not to return nothing (None)

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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    if lamb == None:
        lamb = waveguide.wavelength

    res = ey_1d_x(lamb=lamb, mode=mode, margin=margin, max_solver_iter=max_solver_iter, divisions=divisions, correct_2d=correct_2d)
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
    if no_output:
        return
    return x_vals, field_vals, beta

def plot_ey_1d_y(lamb=None, mode=0, margin=1.5, peak_normalize=False, true_normalize=False, center=False, max_solver_iter=None, divisions=None, correct_2d=False, target_file=None, no_output=False):
    """Evaluates the ey values along the y-axis and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number
    margin : number
        y_vals will be within the range     margin_value * h ± h
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
    correct_2d : bool
        whether or not to produce a more accurate 2d field distribution (though some boundary conditions may not be met)
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot
    no_output : bool
        whether or not to return nothing (None)

    RETURNS
    -------
    (computed values from ey)
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
        - substrate (left)
        - film (middle)
        - cover (right)
    • List of y_vals returned are in terms of h

    """
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    if lamb == None:
        lamb = waveguide.wavelength

    res = ey_1d_y(lamb=lamb, mode=mode, margin=margin, max_solver_iter=max_solver_iter, divisions=divisions, correct_2d=correct_2d)
    if res == None:
        return # return None
    y_vals, field_vals, beta = res
    
    if peak_normalize == True:
        temp = [abs(val) for val in field_vals]
        temp = max(field_vals)
        ec = 0.99 / temp
        field_vals = [val * ec for val in field_vals]
    elif true_normalize == True:
        a = curve_area(y_vals, field_vals)
        ec = 0.9999/a
        field_vals = [val * ec for val in field_vals]
    i = ((int)(2 * max(field_vals)) + 1) / 2
    if center == True:
        y_vals = make_range(-margin, margin, divisions)
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
        
    plt.plot(y_vals, field_vals)

    if target_file != None and type(target_file) == str:
        savefig(target_file)
    if no_output:
        return
    return y_vals, field_vals, beta

def ey_2d_eff(lamb=None, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=None, divisions=None):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    • List of x_vals and y_vals returned are in terms of w and h, respectively

    """
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    res = solve_eff(lamb=lamb, mode=mode, mode_2=mode_2, max_iter=max_solver_iter, ret_trans=True)
    if res == None:
        return None
    beta_prev, beta = res
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

def ey_2d_eff_alt(lamb=None, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=None, divisions=None):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

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

def ey_2d_marcatili(lamb=None, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=None, divisions=None, ignore_beta=False):
    """Evaluates the ey values along both the x and y axes and plots them

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
    margin_x : number
        x_vals will be within the range     margin_value * w ± w
    margin_y : number
        y_vals will be within the range     margin_value * h ± h
    max_solver_iter : number
        max number of iterations to run the solver
    divisions : number
        number of points to evaluate
    ignore_beta : bool
        whether or not to ignore the fact that beta may be less than 0 according to this method

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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

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

def ey_2d_marcatili_ignore_beta(lamb=None, mode=0, mode_2=0, margin_x=10, margin_y=10, max_solver_iter=None, divisions=None):
    """Evaluates the ey values along both the x and y axes and plots them (always ignoring beta)

    PARAMETERS
    ----------
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    • ey_2d_marcatili function such that the fact that beta may be less than 0 is always ignored

    """
    return ey_2d_marcatili(lamb=lamb, mode=mode, mode_2=mode_2, margin_x=margin_x, margin_y=margin_y, max_solver_iter=max_solver_iter, divisions=divisions, ignore_beta=True)

def plot_ey_2d(ey=ey_2d_marcatili, lamb=None, mode=0, mode_2=0, margin_x=2, margin_y=2, max_solver_iter=None, divisions=None, normalize=False, target_file=None, no_output=False):
    """Plots the mode profile looking at a cross-section of the waveguide

    PARAMETERS
    ----------
    ey : func
        chosen function to evaluate ey values
    lamb : number
        wavelength (lambda)
    mode : int
        mode number for transverse propagation in x direction
    mode_2 : int
        mode number for transverse propagation in y direction
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
    no_output : bool
        whether or not to return nothing (None)

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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

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
    if no_output:
        return
    return res

def em_default(lamb=None, mode=0, margin=1.5, max_solver_iter=None, divisions=None):
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
    if max_solver_iter == None:
        max_solver_iter = solvers.default_iter
    if divisions == None:
        divisions = solvers.plot_divisions
    n5, n4, n3, n2, n1 = waveguide.n5, waveguide.n4, waveguide.n3, waveguide.n2, waveguide.n1
    h, w = waveguide.h, waveguide.w
    nf, ns, nc = waveguide.nf, waveguide.ns, waveguide.nc
    if lamb == None:
        lamb = waveguide.wavelength

    start_margin=-margin
    end_margin=margin
    c_light = 3e8
    f = c_light / lamb
    k = 2 * pi / lamb
    V = k * w * sqrt(nf**2 - ns**2)
    a = (ns**2 - nc**2) / (nf**2 - ns**2)
    beta = solve_tm(lamb=lamb, mode=mode, max_iter=max_solver_iter)
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

def plot_em(lamb=None, mode=0, margin=1.5, max_solver_iter=None, divisions=None, target_file=None, no_output=False):
    """Plots the TM mode
    
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
    target_file : string
        target file to save plot to
        if no target file is provided, it does not save the plot to a file and simply returns points to plot
    no_output : bool
        whether or not to return nothing (None)

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

    """
    res = em_default()
    if res == None:
        return 
    plt.plot(res[0], res[1])
    if target_file != None and type(target_file) == str:
        savefig(target_file)
    if no_output:
        return
    return res
