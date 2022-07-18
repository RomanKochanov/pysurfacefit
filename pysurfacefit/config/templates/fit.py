from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'FIT'
    __header__ = 'FITTER OPTIONS'
    __template__ = \
"""
# Weighted / unweighted fit mode.
{weighted_fit}

# Use "rubber" regularization.
{rubber_on}

# Fitting method (trf, lm, basinhopping).
# Available "local" methods:
#    "trf" -> Trust Region Reflective (least squares, bounded).
#    "lm"  -> Levenberg Marquardt (least squares, unbounded).
#    'Nelder-Mead'  -> Gradient-free descent Nelder-Mead method.
#    'Powell'  -> modified Powell algorithm.
#    'CG'  -> conjugate gradient algorithm.
#    'BFGS'  -> BFGS algorithm.
#    'Newton-CG'  -> Newton-CG algorithm.
#    'L-BFGS-B'  -> L-BFGS-B algorithm.
#    'TNC'  -> truncated Newton (TNC) algorithm.
#    'COBYLA'  -> Constrained Optimization BY Linear Approximation.
#    'SLSQP'  -> Sequential Least Squares Programming.
#    'trust-constr'  -> minimize a scalar function subject to constraints.
#    'dogleg'  -> dog-leg trust-region algorithm.
#    'trust-ncg'  -> Newton conjugate gradient trust-region algorithm.
#    'trust-exact'  -> "nearly exact trust-region algorithm"
#    'trust-krylov'  -> "early exact trust-region algorithm"
# Available "global" methods:
#    "basinhopping"  -> local descents from a random starting point (bounded)
#    "anneal" -> Simulated Annealing method
{fitting_method}

# Calculate analytic jacobian (only if ModelSympy is used).
{analytic_jacobian}

# Fit options. Must be given as a list separated by semicolon.
# The valid options for most of the supported fitting methods
# can be found in the corresponding scipy.optimize documentation.  
{fit_options}
"""
    weighted_fit = True
    rubber_on = True
    fitting_method = 'trf'
    analytic_jacobian = True
    fit_options = 'max_nfev=200;'
