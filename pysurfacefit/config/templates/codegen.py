from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'CODEGEN'
    __header__ = 'CODE GENERATOR OPTIONS'
    __template__ = \
"""
# Create Fortran code.
{create_fortran}

# Compare original Python code with the generated Fortran.
{compare_fortran}

# Fortran compiler executable
{compiler_fortran}

# Grid specifications to calculate on.
# Format of the specification must be as follows: 
# X=XMIN:XSTEP:XMAX; Y=YMIN:YSTEP:YMAX; ...
{gridspec}
"""
    create_fortran = False
    compare_fortran = False
    compiler_fortran = 'ifort'
