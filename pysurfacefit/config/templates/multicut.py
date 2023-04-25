from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'MULTICUT'
    __header__ = 'MULTICUT OPTIONS'
    __template__ = \
"""
# Plot multiple 1D cuts of the datapoints and functions. 
# To be used when spotting the outliers, and imaging the behaviour
# of the fitted model.

# If there are too many cuts, setting this parameter to False can 
# save the screen space.
{show_legend}

# Draw lines for the data points. Setting to False will make the space
# less busy when there are too many cuts.
{show_lines}

# Plotting model flag.
{plot_model}

# If set to True, model (main and aux.) will be calculated on a full
# grid. If set to False, they will be calculated only where the data 
# are available.
{fullgrid_calc}

# Auxiliary models to compare with. Should be splitted with a semicolon.
{compare_with_models}

# Text field to display with data points. Should be an existing field 
# in all datafiles.
{annotate_field}

# Unfixed argument along the x axis.
{argument}

# Gridspec definind the bounds and grids for filtering data and
# calculating models. 
{gridspec}
"""
    show_legend = True
    show_lines = False
    plot_model = True
    fullgrid_calc = False
