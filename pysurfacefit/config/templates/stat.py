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

# Outlier statistics global flag.
# If set to True, statistics are calculated on the full Jacobian.
# If set to False, the Jacobian is split into datagroups and all stats are 
# calculated separately for each datagroup.
{outlier_stats_global}

# Type of statistics: cook, dffits, leverage ,student.
{outlier_stats_type}

# Output generated symbolic function, for debug purposes only.
{output_symbolic_func}
"""
    stat_file = 'stat.out'
    outlier_stats_flag = True
    outlier_stats_global = True
    outlier_stats_type = 'cook'
