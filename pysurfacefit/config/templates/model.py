from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'MODEL'
    __header__ = 'FITTING MODEL SPECIFICATIONS'
    __template__ = \
"""
# Model package name.
{model}

# Model arguments. Argument names must be in the same order
# as the data column names from the DATA section.
# E.g.: X;Y;Z
{arguments}
"""
    model = 'fitmodel'
