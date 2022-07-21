from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'CALC'
    __header__ = 'CALCULATION OPTIONS'
    __template__ = \
"""
# Output file to output calculated model. 
# Column names are defined in DATA section.
{output_file}

# Grid specifications to calculate on.
# Format of the specification must be as follows: 
# X=XMIN:XSTEP:XMAX; Y=YMIN:YSTEP:YMAX; ...
{gridspec}
"""
    output_file = 'calc.csv'