from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'MODEL'
    __header__ = 'FITTING MODEL SPECIFICATIONS'
    __template__ = \
"""
# Model package name to generate.
{model}

# Model arguments. Argument names must be in the same order
# as the data column names from the DATA section.
# E.g.: X;Y;Z
{arguments}

# Name of the parameter file to use in fitting.
{parfile}
"""
    model = 'fitmodel'
    parfile = 'fit.par'
