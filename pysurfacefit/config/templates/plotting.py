from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'PLOTTING'
    __header__ = 'PLOTTING OPTIONS'
    __template__ = \
"""
# ATTENTION: names of the plot coordinates correspond to the 
# model argumen names given in the MODEL section.

# Type of the plot.
# Available plot modes:
#   => "residuals": 
#           plot unweighted fit residuals.
#   => "sections": 
#           plot sections/cuts of the fitted model vs datapoints.
#   => "history":
#           plot fitting history for each data group.
#   => "pareto":
#           plot Pareto frontier for the Leverage/Cook outlier pairs.
{plot_mode}

# MIN_WHT_RES; MAX_WHT_RES; MIN_UNWHT_RES; MAX_UNWHT_RES; WHT_SSE; UNWHT_SSE; WHT_SD; UNWHT_SD
{fit_history_stat}
{fit_history_logscale} 

# Plot coordinate grid specifications.
# Format of the specification must be as follows: 
# X=XVAL; Y=YMIN:YSTEP:YMAX; ...
# In case of the fixed coordinate, variable would only have a value, e.g. X=XVAL.
# In case of unfixed coordinate, there should be a full grid, e.g. Y=YMIN:YSTEP:YMAX
# N.B.: 1) order of the unfixed coords affects the order of the plot axes.
#       2) order of the binding MUST correspond to the argument of the model's 
#          __func__ method.
{gridspec}

# Model components to plot.
{model_components}

# Calculate model's components.
{calculate_components}

# Plot outlier statistics in color. If False, each datagroup has its own color.
{plot_outlier_stats}

# Plot weighted residuals.
{resids_weighted}

# X axes to plot the resuduals versus. 
# If empty, defaults to [DATA][output_column].
{resids_x_axes}

# Scatter settings (2D, 3D case)
{scatter_opacity}
{marker_size}
{resize_by_weights}

# Surface settings (3D case)
{surface_opacity}
"""
    plot_mode = 'residuals'
    fit_history_stat = 'WHT_SSE'
    fit_history_logscale = True
    scatter_opacity = 1
    resize_by_weights = False
    surface_opacity = 0.4
    calculate_components = False
    resids_weighted = False
    plot_outlier_stats = False
    marker_size = 20
