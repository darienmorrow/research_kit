import numpy as np
import WrightTools as wt
from scipy.optimize import least_squares as ls 
import time 


def simple_fit(objective_func, guess_params, xdata, ydata, verbose=True, bounds=(-np.inf, np.inf)):
    """
    Wrapper for SciPy's scipy.optimize.least_squares

    Parameters
    ----------
    objective_func : function
        function to fit to with syntax of f(x, p1, p2, ...,) or f(x, params)
    guess_params : array_like with shape (n,)
        Initial guess on independent variables. 
    xdata : array_like with shape (m,)
        Independent data.
    ydata : array_like with shape (m,)
        Dependent data.
    verbose : bool, optional
        Toggle printing of fitting time, cost, and parameters. The default is True.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. 
        Defaults to no bounds. 
        Each array must match the size of guess_params or be a scalar, in the latter case a bound will be the same for all variables. 
        Use np.inf with an appropriate sign to disable bounds on all or some variables.
        The default is (-np.inf, np.inf).

    Returns
    -------
    params : array_like with shape (n,)

    """
    def cost_func(params):
        # simple, linear cost function. 
        # One can make more nuanced cost functions if needed.
        cost = ydata - objective_func(xdata, params)
        return cost
    t0 = time.time()
    out = ls(cost_func, guess_params, bounds=bounds)
    t1 = time.time()
    if verbose:
        print('Fit done in {0} s with cost of {1}.'.format(str(round(t1-t0,2)), str(out.cost)))
        print(out.x)
    return out.x


def gauss(t, t0, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))


def exp(t, t1, A1, B, t0):
    # applies a heaviside step function
    zero = t0
    out = np.zeros(t.size)
    out[t >= zero] = A1 * np.exp(-t[t >= zero] / t1)
    out[t == zero] *= 0.5
    out += B
    return out


def biexp(t, t1, A1, t2, A2, B, t0):
    # applies a heaviside step function
    zero = t0
    out = np.zeros(t.size)
    out[t >= zero] = A1 * np.exp(-t[t >= zero] / t1) + A2 * np.exp(-t[t >= zero] / t2)
    out[t == zero] *= 0.5
    out += B
    return out


def triexp(t, t1, A1, t2, A2, t3, A3, B, t0):
    # applies a heaviside step function
    zero = t0
    out = np.zeros(t.size)
    out[t >= zero] = (
        A1 * np.exp(-t[t >= zero] / t1)
        + A2 * np.exp(-t[t >= zero] / t2)
        + A3 * np.exp(-t[t >= zero] / t3)
    )
    out[t == zero] *= 0.5
    out += B
    return out


def exp_fit_func(p, x):
    # oversample grid. Then convolve, then interpolate convolution back to original grid
    t = np.linspace(
        x.min(), x.max(), 1024 * 4
    )  # need uniformly spaced grid to do convolution
    t1, A1, B, t0, fwhm = p
    t1, A1, B, fwhm = np.abs(t1), np.abs(A1), np.abs(B), np.abs(fwhm)
    IRF = gauss(t - t.mean(), t0, fwhm)
    IRF /= IRF.sum()  # need area normalized function
    decay = exp(t, t1, A1, B, t0)
    model = np.convolve(decay, IRF, mode="same")
    out = np.interp(x, t, model)
    return out


def biexp_fit_func(p, x):
    # oversample grid. Then convolve, then interpolate convolution back to original grid
    t = np.linspace(
        x.min(), x.max(), 1024 * 4
    )  # need uniformly spaced grid to do convolution
    t1, A1, t2, A2, B, t0, fwhm = p
    t1, A1, t2, A2, B, fwhm = (
        np.abs(t1),
        np.abs(A1),
        np.abs(t2),
        np.abs(A2),
        np.abs(B),
        np.abs(fwhm),
    )
    IRF = gauss(t - t.mean(), t0, fwhm)
    IRF /= IRF.sum()  # need area normalized function
    decay = biexp(t, t1, A1, t2, A2, B, t0)
    model = np.convolve(decay, IRF, mode="same")
    out = np.interp(x, t, model)
    return out


def triexp_fit_func(p, x):
    # oversample grid. Then convolve, then interpolate convolution back to original grid
    t = np.linspace(
        x.min(), x.max(), 1024 * 4
    )  # need uniformly spaced grid to do convolution
    t1, A1, t2, A2, t3, A3, B, t0, fwhm = p
    t1, A1, t2, A2, t3, A3, B, fwhm = (
        np.abs(t1),
        np.abs(A1),
        np.abs(t2),
        np.abs(A2),
        np.abs(t3),
        np.abs(A3),
        np.abs(B),
        np.abs(fwhm),
    )
    IRF = gauss(t - t.mean(), t0, fwhm)
    IRF /= IRF.sum()  # need area normalized function
    decay = triexp(t, t1, A1, t2, A2, t3, A3, B, t0)
    model = np.convolve(decay, IRF, mode="same")
    out = np.interp(x, t, model)
    return out


def exp_param_guess(x, y):
    t0 = x[np.argmax(y)]
    B = np.mean(y[-100:])
    fwhm = (x.min() + t0) / 2
    A1 = y.max() / 2
    t1 = np.mean(x)
    p = [t1, A1, B, t0, fwhm]
    return p


def biexp_param_guess(x, y):
    t0 = x[np.argmax(y)]
    B = np.mean(y[-100:])
    fwhm = (x.min() + t0) / 2
    A1 = y.max() / 2
    t1 = np.mean(x)
    A2 = A1 / 5
    t2 = t1 * 10
    p = [t1, A1, t2, A2, B, t0, fwhm]
    return p


def triexp_param_guess(x, y):
    t0 = x[np.argmax(y)]
    B = np.mean(y[-100:])
    fwhm = (x.min() + t0) / 2
    A1 = y.max() / 2
    t1 = np.mean(x)
    A2 = A1 / 5
    t2 = t1 * 10
    A3 = A1 / 10
    t3 = t1 * 100
    p = [t1, A1, t2, A2, t3, A3, B, t0, fwhm]
    return p


def sqrt_fit(p0, x, y, func):
    sqrty = np.sqrt(y)

    def sqrtfunc(p, x):
        return np.sqrt(func(p, x))

    pfit, perr = wt.kit.leastsqfitter(p0, x, sqrty, sqrtfunc)
    return pfit, perr


def exp_fit(x, y):
    p0 = exp_param_guess(x, y)
    pfit, perr = sqrt_fit(p0, x, y, exp_fit_func)
    ymodel = exp_fit_func(pfit, x)
    return pfit, perr, ymodel


def biexp_fit(x, y):
    p0 = biexp_param_guess(x, y)
    pfit, perr = sqrt_fit(p0, x, y, biexp_fit_func)
    ymodel = biexp_fit_func(pfit, x)
    return pfit, perr, ymodel


def triexp_fit(x, y):
    p0 = triexp_param_guess(x, y)
    pfit, perr = sqrt_fit(p0, x, y, triexp_fit_func)
    ymodel = triexp_fit_func(pfit, x)
    return pfit, perr, ymodel
