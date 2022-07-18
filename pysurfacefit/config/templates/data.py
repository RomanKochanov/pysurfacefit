from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'DATA'
    __header__ = 'FITTING DATA SPECIFICATION'
    __template__ = \
"""
# CSV file containing fitting data.
{datafile}

# Data specification file of the following format:
# alias   path       wht_fun    filter   include
# ------- ---------  ---------  -------  -------
# test1   test1.csv  wht1       filter1  0       
# test2   test2.csv  wht2       1
# test3   test3.csv  wht3       1
{dataspec}

# Names and units of the input data columns. 
# Input names must be separated by semicolon 
# and should not contain dots.
# E.g.: X;Y;Z for names, X:X_UNIT;Y:Y_UNIT for units
{input_columns}
{input_units}

# Name of the output data column. See explanation above.
{output_column}
{output_units}

# Weight function.
{wht_fun}

# Global data cutoff
{global_cutoff_max}
{global_cutoff_min}
"""
    wht_fun = 'lambda *args: 1'
    dataspec = 'dataspec.txt'
