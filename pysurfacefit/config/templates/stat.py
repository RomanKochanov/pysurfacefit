from ..base import ConfigSection

class Default(ConfigSection):
    __name__ = 'STAT'
    __header__ = 'STATISTICS CALCULATION OPTIONS'
    __template__ = \
"""
# Fit statistics file.
{stat_file}

# Calculate outlier statistics.
{outlier_stats_flag}

# Type of statistics: cook, dffits, leverage ,student.
{outlier_stats_type}

# Output generated symbolic function, for debug purposes only.
{output_symbolic_func}
"""
    stat_file = 'stat.out'
    outlier_stats_flag = True
    outlier_stats_type = 'cook'
