from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'DATA'
    __header__ = 'FITTING DATA SPECIFICATION'
    __template__ = \
"""
# CSV file containing fitting data.
{datafile}

# Define rules to split the datafile to dataspec.
# Should be ignored if dataspec is empty.
{split_column}
{split_values}
{split_weights}

# Resulting data specification split-file.
# N.B.: if empty, then the raw datafile is used.
# Data specification *.txt file has the following format:
# alias   path       wht_mul   type    include
# ------- ---------  --------  ------  -------
# test1   test1.csv  1         0       1       
# test2   test2.csv  1         0       1
# test3   test3.csv  1         0       1
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
    wht_fun = 'lambda v: 1'
    dataspec = 'dataspec.txt'
