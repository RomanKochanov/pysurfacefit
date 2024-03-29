import os
import re
import sys
import json
import shutil

import numpy as np
import jeanny3 as j

from functools import reduce
from itertools import cycle

from .data import FitGroups, FitPoints, PenaltyPointsUp, PenaltyPointsDown
from .grids import Grid, List
from .fitter import Fitter
from .formats.dill import serialize_fit, deserialize_fit 
from .formats.dill import serialize_model, deserialize_model
from .models.sympy import create_pes_interface, create_pots_interface

import matplotlib 
matplotlib.use('TkAgg') # backend should be set before pylab and pyplotl are imported!
# 'GTK', 'GTKAgg', 'GTKCairo', 'GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 
# 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 
# 'WXCairo', 'agg', 'cairo', 'gdk', 'pdf', 'pgf', 'ps', 'svg', 'template'
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

MARKERS = ['o','v','^','<','>','s','p','P','*','+','x','D','d']

#COLORS = ['blue', 'green', 'red', 'black', 'cyan', 
#'magenta', 'indigo', 'chartreuse', 'yellow', 'chocolate', 'aqua', 'aquamarine', 'azure', 'beige', 'brown', 
#'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'goldenrod', 
#'grey', 'ivory', 'coral', 'khaki', 'lavender', 'lightblue', 'lightgreen', 'lime', 
#'maroon', 'navy', 'olive', 'orange', 'orangered', 'orchid', 'pink', 'plum', 
#'purple', 'salmon', 'sienna', 'silver', 'tan', 'teal', 'tomato', 'turquoise', 
#'violet', 'wheat', 'yellowgreen']

COLORS = [
'blue', 'green', 
#'red', 
'black', 'cyan', 
#'magenta', 'indigo', 'chartreuse', 'yellow', 'chocolate', 'aqua', 'aquamarine', 'azure', 'beige', 'brown',
'magenta', 'indigo', 'chartreuse', 'yellow', 'aqua', 'aquamarine', 'brown', 
#'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'goldenrod', 
'darkblue', 'darkgreen', 'fuchsia', 'gold', 'goldenrod', 
'grey', 
#'ivory', 
'coral', 'khaki', 'lavender', 'lightblue', 'lightgreen', 'lime', 
'maroon', 'navy', 'olive', 'orange', 'orangered', 'orchid', 'pink', 'plum', 
'purple', 'salmon', 'sienna', 'silver', 'tan', 'teal', 'tomato', 'turquoise', 
'violet', 'wheat', 'yellowgreen']

COLORS_SUPPORT = [
    'red', 
    'magenta', 'indigo', 'chartreuse', 'yellow', 'chocolate', 'aqua', 'aquamarine', 'azure', 'beige', 'brown',
    'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'goldenrod', 
]

# type conversions minding the None value
to_type = lambda val,typ: typ(val) if val not in (None,'') else None
to_int = lambda val: to_type(val,int)
to_float = lambda val: to_type(val,float)
#to_bool = lambda val: val if val.lower() in ("yes", "true", "t", "1")
def to_bool(val):
    if type(val) is bool:
        return val
    elif val is None:
        return False
    elif type(val) in {int,float}:
        return bool(val)
    elif type(val) is str:
        return val.lower() in ("yes", "true", "t", "1")
    else:
        raise Exception('unknown type')

######################
#### STARTPROJECT ####
######################

def create_dir(dirname):
    """ Create a directory or complain if existing"""
    if os.path.isdir(dirname):
        print('Error: "%s" directory already exists'%dirname)
        sys.exit()
    os.mkdir(dirname)

def create_file(filename,content):
    if os.path.isfile(filename):
        print('Error: "%s" file already exists'%filename)
        sys.exit()
    with open(filename,'w') as f:
        f.write(content)

def create_initfile(subdir,content):
    """ Create __init__.py file for a Python package"""
    with open(os.path.join(subdir,'__init__.py'),'w') as f:
        f.write(content)

def create_dataspec(CONFIG):
    """ Create a sample Data Specification file. """
    project_dir = CONFIG['GENERAL']['project']
    filename = CONFIG['DATA']['dataspec']
    datafile = CONFIG['DATA']['datafile']
    content = """//HEADER
A alias str
P path str
W wht_mul float
T type int
I include int

//DATA
A_________________P_______________________W__________T__________I______
# type: 0/None -> data, -1 -> penalty down, 1 -> penalty up
"""
    if datafile:
        filestem,_ = os.path.splitext(datafile)
        content += """
{filestem}            {datafile}            1.0         0         1
""".format(filestem=filestem,datafile=datafile)
    #with open(os.path.join(project_dir,filename),'w') as f:
    #    f.write(content)
        
    col = j.Collection()
    col.update([{
        'alias': filestem,
        'path': datafile,
        'wht_mul': 1.0,
        'type': 0,
        'include': 1,
    }])
    col.order = ['alias','path','wht_mul','type','include']
    col.export_fixcol(os.path.join(project_dir,filename))
        
    return filename

def startproject(CONFIG):
    project_dir = CONFIG['GENERAL']['project']
    datafile = CONFIG['DATA']['datafile']
    create_dir(project_dir)
    if datafile:
        shutil.copy(datafile,os.path.join(project_dir,datafile)) 
    config_name = 'config.ini'
    config_path = os.path.join(project_dir,config_name)
    CONFIG.save_ini(config_path)
    dataspec_file = create_dataspec(CONFIG)
    print('Created new project: %s'%project_dir)
    print('New config file has been added: %s'%config_path)
    print('Sample data specification file has been added: %s'%dataspec_file)

##############
#### INIT ####
##############

def parse_input_columns(CONFIG):
    colnames = [c.strip() for c in CONFIG['DATA']['input_columns'].split(';')]
    if colnames[-1]=='': colnames = colnames[:-1]
    return colnames

def parse_semicolon_list(CONFIG,section,parameter,typ=str):
    vals = [typ(c.strip()) for c in CONFIG[section][parameter].split(';')]
    if vals[-1]=='': vals = vals[:-1]
    return vals

def split(CONFIG):
    """ Split datafile according to the split scheme. """
    
    # Get the config options.
    datafile = CONFIG['DATA']['datafile'].strip()
    dataspec = CONFIG['DATA']['dataspec'].strip()
    if not dataspec: dataspec = 'dataspec.txt'
    split_column = CONFIG['DATA']['split_column'].strip()
    if not split_column: 
        split_column = CONFIG['DATA']['output_column'].strip()
    split_values = parse_semicolon_list(CONFIG,'DATA','split_values',typ=float)
    split_weights = parse_semicolon_list(CONFIG,'DATA','split_weights',typ=float)
    
    print('Splitting data by %s'%split_column)
    
    col_data = j.import_csv(datafile)
    
    filename = os.path.basename(datafile)
    filepath = os.path.dirname(datafile)
    filestem,_ = os.path.splitext(filename)
        
    if not split_values:
        cols = [col_data]
        whts = [split_weights[0]] if split_weights else [1.0]
    else:
        cols = col_data.split(split_column,split_values)
        whts = split_weights + [1.0 for _ in range(len(split_values)-len(split_weights))]
        whts += [whts[-1]]
        split_values += [None]
        assert len(whts)==len(split_values)==len(split_values)

    col_dataspec = j.Collection()

    _,dataspec_ext = os.path.splitext(dataspec)
    if dataspec_ext.lower()=='.txt':
        dataspec_save_method = col_dataspec.export_fixcol
    elif dataspec_ext.lower()=='.csv':
        dataspec_save_method = col_dataspec.export_csv
    else:
        print('ERROR: unknown dataspec export format "%s"'%dataspec_ext)
        sys.exit()
    
    print('Saving dataspec to ',dataspec)

    if len(cols)==1:
        col_dataspec.update({
            'alias':filestem,
            'path':datafile,   
            'wht_mul':whts[0],    
            'type':None,
            'include':1,
        })
    else:
        for i in range(len(split_values)):
            
            col_ = cols[i]
            
            npnts = len(col_.ids())
            if npnts==0:
                continue
            
            if i==0:
                filestem_ = filestem+'__%g-%g'%(
                    min(col_.getcol(split_column)),
                    split_values[i],
                )
            elif i==len(split_values)-1:
                filestem_ = filestem+'__%g-%g'%(
                    split_values[i-1],
                    max(col_.getcol(split_column)),
                )
            else:
                filestem_ = filestem+'__%g-%g'%(
                    split_values[i-1],
                    split_values[i],
                )
                
            datafile_ = os.path.join(filepath,filestem_+'.csv') 
            
            col_dataspec.update({
                'alias': filestem_,
                'path': datafile_,   
                'wht_mul': whts[i],    
                'type': None,
                'npnts': npnts,
                'include': 1,
            })
            
            col_.export_csv(datafile_)
    
    col_dataspec.order = ['alias','path','wht_mul','type','npnts','include']
    dataspec_save_method(dataspec)
    print('...done')

def parse_dataspec(CONFIG):
    filename = CONFIG['DATA']['dataspec'].strip()
    if filename:
        col_dataspec = j.import_fixcol(filename)
    else:
        col_dataspec = j.Collection()
    return col_dataspec

def read_fitgroups(CONFIG,verbose=0,exclude=True):
    """ Read data fitgroups based on the config file information. """
    
    # Create unified collection for plotting.
    col_data = j.Collection() # unified collection for plotting

    # Get the config options.
    datafile = CONFIG['DATA']['datafile'].strip()
    dataspec = CONFIG['DATA']['dataspec'].strip()
    
    exclude_file = CONFIG['DATA']['exclude']
    if exclude_file is not None:
        exclude_file = exclude_file.strip()
    
    # => input and output names
    INPUTS = parse_input_columns(CONFIG)
    OUTPUT = CONFIG['DATA']['output_column']
    
    # => cutoff min and max values
    GLOBAL_CUTOFF_MIN = to_float( CONFIG['DATA']['global_cutoff_min'] )
    GLOBAL_CUTOFF_MAX = to_float( CONFIG['DATA']['global_cutoff_max'] )
    if GLOBAL_CUTOFF_MIN is None: GLOBAL_CUTOFF_MIN=-np.inf
    if GLOBAL_CUTOFF_MAX is None: GLOBAL_CUTOFF_MAX=np.inf
    
    # => weight functin for data
    wht_fun = eval( CONFIG['DATA']['wht_fun'] )
    
    # Import exclusion collection and make a lookup hash
    if exclude_file and exclude:
        col_exclude = j.import_csv(exclude_file)
        #get_exclude_keys = lambda v: tuple([v[k] for k in INPUTS+[OUTPUT]])
        get_exclude_keys = lambda v: tuple([v[k] for k in INPUTS])
        exclude_hash = col_exclude.group(get_exclude_keys)
    
    if not datafile and not dataspec:
        print('ERROR: either datafile of dataspec should be non-empty.')
        sys.exit()
    
    col_dataspec = j.Collection()
    # If dataspec is empty, treat datafile as a single data source:
    if datafile and not dataspec:
        if verbose>0: print('Reading from datafile',datafile)
        col_dataspec.update({
            'alias':'default','path':datafile,
            'wht_mul':1.0,'type':None,'include':1})
    
    # Second, get fitgroups from the data in the csv file.
    col = parse_dataspec(CONFIG)
    col_dataspec.update(col.getitems())
    
    order = []
    
    col_summary = j.Collection()
    col_summary.order = ['group','N_points','excluded','type']
    
    GROUPS_DATA = []
    for item in col_dataspec.getitems():
        item_summary = {}
        # Get fitgroup properties
        fitgroup_path = item['path'].strip()
        fitgroup_name = item['alias'].strip()
        fitgroup_type = item['type']
        fitgroup_wht_mul = item['wht_mul']
        item_summary['group'] = fitgroup_name
        # Get points from the CSV file
        col = j.import_csv(fitgroup_path)
        order = col.order
        # Apply data cutoff
        col = col.subset(col.ids(lambda v: GLOBAL_CUTOFF_MIN<=v[OUTPUT]<=GLOBAL_CUTOFF_MAX)); n = len(col.ids())
        # Apply data exclusion
        if exclude_file and exclude: 
            col = col.subset(col.ids(lambda v: get_exclude_keys(v) not in exclude_hash))
        item_summary['N_points'] = n
        item_summary['excluded'] = n-len(col.ids())
        # Get the data columns 
        input_columns = col.getcols(INPUTS)
        output_column = col.getcol(OUTPUT)
        # Calculate data weights 
        col.assign('fit_weight',wht_fun)
        data_whts = col.getcol('fit_weight')
        # Print collection contents
        if verbose>=2:
            print('\nGROUP NAME',fitgroup_name)
            col.tabulate(['N',*(INPUTS+[OUTPUT]),'fit_weight'])
        # Prepare fitting grid (list)
        grid = List(*reduce(lambda x,y:x+y,zip(INPUTS,input_columns))) # list-based grid
        # Prepare fitgroup for points
        if type(fitgroup_type) is str:
            fitgroup_type = fitgroup_type.strip()
            if fitgroup_type: fitgroup_type = int(fitgroup_type)
        if fitgroup_type in ['', None, 0]:
            FitPointsType = FitPoints
        elif fitgroup_type==-1:
            FitPointsType = PenaltyPointsDown
        elif fitgroup_type==1:
            FitPointsType = PenaltyPointsUp
        else:
            print('ERROR: unknown dataspec type: %s'%str(fitgroup_type))
            sys.exit()
        item_summary['type'] = FitPointsType.__name__
        GRP_DATA = FitPoints(
            name=fitgroup_name,input=grid,output=output_column,
            whtmul=fitgroup_wht_mul,whts=data_whts,active=True
        )
        # Add to the plotting collection
        col_data.update(col.getitems())
        # Append to groups of data
        if col.ids():
            GROUPS_DATA.append(GRP_DATA)
        else:
            if verbose>=1: 
                print(' // skipping empty group(%s)'%fitgroup_name)
        
        # Update summary
        col_summary.update(item_summary)
        
    if verbose>0:
        print('')
        print('==============================')
        print('===== FITGROUPS SUMMARY ======')
        print('==============================')
        print('')
        col_summary.tabulate()
               
    col_data.order = order
               
    fitgroups = FitGroups(GROUPS_DATA)
    return fitgroups,col_data

def get_data_units(CONFIG):
    """ Get units of the data points. Used in codegen and model creation. """
    
    input_unitspec = CONFIG['DATA']['input_units']
    output_unitspec = CONFIG['DATA']['output_units']
    
    data_inputs = CONFIG['DATA']['input_columns']
    model_args = CONFIG['MODEL']['arguments']
    input_mapping = {}
    for arg,colname in zip(model_args.split(';'),data_inputs.split(';')):
        arg = arg.strip(); colname = colname.strip()
        input_mapping[colname] = arg
    
    input_units = {}
    if input_unitspec:
        for pair in input_unitspec.split(';'):
            pair = pair.strip(); 
            if not pair: continue
            argname,units = pair.split(':')
            argname = argname.strip()
            units = units.strip()
            input_units[input_mapping[argname]] = units
    
    if output_unitspec:
        argname,units = output_unitspec.split(':')
        output_units = output_unitspec
    else:
        output_units = 'None'
    
    return {
        'input': input_units,
        'output': output_units
    }

def create_model_package(CONFIG):
    module_name = CONFIG['MODEL']['model']
    #model_name = 'Model'
    model_name = module_name
        
    template = """import os

from pysurfacefit.models.sympy import ModelSympy
from pysurfacefit.fitpars import Par, Parameters
    
class {modelname}(ModelSympy):
    \"\"\"
    Test model for fitting data with just one parameter.
    \"\"\"
    def __init__(self,calc_switch='numbified'):
        self.__check_symbolic__ = False
        self.__calc_switch__ = calc_switch # symbolic, lambdified, numbified
        self.__components__ = {}
        
        # Initialize empty parameters.
        self.__params__ = Parameters()
               
        # Add constant parameter.
        self.__params__.append(group='constant', pars=[
            Par(name='constant', value=0.0, flag=True),
        ])
        
    def __units__(self):
        return {dataunits}
        
    def __func__(self,params,{arguments}):
                
        # constant
        constant = params['constant'].get_value()
                
        # final result
        res = constant
        
        # components
        self.__components__['constant'] = constant
                    
        return res
        
model = {modelname}()

if os.path.exists('{modelname}.csv'):
    model.load_params('{modelname}.csv')
else:
    model.save_params('{modelname}.csv')
"""
     
    content = template.\
        replace('{modelname}',model_name).\
        replace('{dataunits}',json.dumps(get_data_units(CONFIG))).\
        replace('{arguments}',','.join(get_arguments(CONFIG)))
     
    # Create a "one-file" model package.
    create_file(module_name+'.py',content)

def init(CONFIG):
    read_fitgroups(CONFIG)
    create_model_package(CONFIG)

##############
#### FIT #####
##############

def parse_fit_options(CONFIG):
    """ Parse fit options from the FIT config section.
        Options shoudl be given in the following format:
        opt1=val1;opt2=val2;...
        Automatic type inference is applied."""    
    fit_options_line = CONFIG['FIT']['fit_options']
    pairs = fit_options_line.split(';')
    fit_options = {}
    for pair in pairs:
        if not pair: continue
        name,value = pair.split('=')
        name = name.strip()
        value = value.strip()
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
        fit_options[name] = value
    return fit_options

def import_module_by_path(module_name,module_path):
    # https://www.geeksforgeeks.org/how-to-import-a-python-module-given-the-full-path/
    # https://stackoverflow.com/questions/67631/how-do-i-import-a-module-given-the-full-path
    import importlib.util
    spec=importlib.util.spec_from_file_location(module_name,module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module # dill/pickle doesn't work without this line
    spec.loader.exec_module(module)
    return module

def load_module_(module_name):
    """ Load model from its module (no config).
        Model should read its parameters automatically! """
    module_path = os.path.join('./',module_name+'.py') # relative path
    module = import_module_by_path(module_name,module_path)
    return module

def load_model_(module_name):
    """ Load model from its module (no config).
        Model should read its parameters automatically! """
    module = load_module_(module_name)
    model = getattr(module,'model')
    return model

def load_module_rec(module_name,verbose=0):
    """ Recursively load module.
        Fixes errors with simple load_model_ when there are nested imports
        of the local modules inside the model script. """
    while True:
        try:
            if verbose>0: print('Trying to import',module_name)
            module = load_module_(module_name)
            if verbose>0: print('  success!')
            flag_success = True
        except ModuleNotFoundError as e:
            missing_module_name = e.name
            if verbose>0: print('  failure... Need to import %s first.'%missing_module_name)
            if missing_module_name==module_name:
                print('ERROR: cannot find module "%s"'%missing_module_name)
                sys.exit()
            load_module_rec(missing_module_name,verbose)
            flag_success = False
        if flag_success:
            if verbose>0: print('breaking for',module_name)
            break
    return module
    
def load_model_rec(module_name,verbose=0):
    module = load_module_rec(module_name,verbose)
    model = getattr(module,'model')
    return model

def load_model(CONFIG,verbose=0):
    """ Load model from its module.
        Model should read its parameters automatically! """
    module_name = CONFIG['MODEL']['model']
    # Old non-recursive import, fails when there is a nested local imports 
    # inside the model script:
    #model = load_model_(module_name) 
    # New recursive import:
    #verbose = 1 # debug
    if verbose>0:
        print('')
        print('================================')
        print('Starting import for',module_name)
        print('================================')
    model = load_model_rec(module_name,verbose)
    return model

def fit(CONFIG):
    """ Start fitting procedure. """
    
    # Get options from the config file.
    weighted_fit = to_bool(CONFIG['FIT']['weighted_fit'])
    rubber_on = to_bool(CONFIG['FIT']['rubber_on'])
    fitting_method = CONFIG['FIT']['fitting_method']
    analytic_jacobian = to_bool(CONFIG['FIT']['analytic_jacobian'])
    stat_file = CONFIG['STAT']['stat_file']
    fitopts = parse_fit_options(CONFIG)
    model_name = CONFIG['MODEL']['model']
    model = load_model(CONFIG)
    fitgroups,_ = read_fitgroups(CONFIG,verbose=1)
    
    # Create Fitter object from scratch.
    f = Fitter(model=model,fitgroups=fitgroups,
            weighted=weighted_fit,rubber_on=rubber_on,
            method=fitting_method,jac=analytic_jacobian,
            **fitopts);
            
    # Start fit.
    f.fit()
    
    # Save model.
    f.model_final.save_params(model_name+'.csv')
    
    # Get Jacobian, if any.
    f.__fitgroups__.__jac__ = getattr(f.__result__,'jac',None)
    
    # Save statistics.
    with open(stat_file,'w') as stat_stream:
        calculate_and_save_statistics(
            CONFIG,f.model_final,f.__fitgroups__,
            stream=stat_stream,active_only=True)
    
##############
#### STAT ####
##############

def calculate_fitpar_statistics(CONFIG,fitgroups,save=True,stream=sys.stdout,active_only=True):
    
    model_name = CONFIG['MODEL']['model']
    
#    # Output covariance matrix
#    print('\n10-ORDERS OF COVARIANCE ("!" MEANS NEGATIVE VALUE)')
#    cc = j.Collection(); 
#    for i,row in enumerate(f.__covb__):
#        #item = {str(j):e for j,e in enumerate(row)}
#        item = {str(j):(str(int(np.log10(abs(e))))+\
#            ('' if e>0 else '!')) if j>=i else None for j,e in enumerate(row)} # show only orders of covariance
#        item['#'] = i
#        cc.update(item)
#    order = ['#']+[str(i) for i in range(f.__covb__.shape[0])]
#    cc.tabulate(order)

def calculate_group_statistics(CONFIG,fitgroups,save=True,stream=sys.stdout,active_only=True):
    
    model_name = CONFIG['MODEL']['model']

#https://mlink.in/qa/?qa=769607/
#https://code.activestate.com/recipes/578287-multidimensional-pareto-front/
def is_pareto_efficient(costs,maximize=True,return_mask=True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        if maximize:
            nondominated_point_mask = np.any(costs>costs[next_point_index], axis=1)
        else:
            nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient
        
def calculate_residual_statistics(CONFIG,fitgroups,save=True,stream=sys.stdout,active_only=True):
    
    model_name = CONFIG['MODEL']['model']

    output_column = CONFIG['DATA']['output_column']
    
    output = fitgroups.get_output(active_only=active_only)
    input_matrix = fitgroups.get_inputs(active_only=active_only)
    full_weights = fitgroups.get_weights(active_only=active_only)

    #model_calc = fitgroups.collect('__calc_vals__',active_only=active_only)
    #weighted_residuals = fitgroups.collect('__weighted_residuals__',active_only=active_only)
    #unweighted_residuals = fitgroups.collect('__unweighted_residuals__',active_only=active_only)
    model_calc = fitgroups.__calc_vals__
    weighted_residuals = fitgroups.__weighted_resids__
    unweighted_residuals = fitgroups.__unweighted_resids__
    
    length = fitgroups.get_length()
    input_names = fitgroups.__vars__
    
    outlier_stats_flag = to_bool( CONFIG['STAT']['outlier_stats_flag'] )
    outlier_stats_global = to_bool( CONFIG['STAT']['outlier_stats_global'] )
    
    #print('===============calculate_residual_statistics=====================')
    #print('outlier_stats_flag>>>',outlier_stats_flag,type(outlier_stats_flag))
    #print('outlier_stats_global>>>',outlier_stats_global,type(outlier_stats_global))
    
    outlier_stats = []
    if outlier_stats_flag and fitgroups.__jac__ is not None: 
        
        #print('CALCULATING OUTLIER STATS:')
        
        if outlier_stats_global:
            STATS = fitgroups.calculate_outlier_statistics()
            for stat in STATS:
                fitgroups.split(STATS[stat],stat,active_only=True)
        else:
            fitgroups.split_global_jacobian(fitgroups.__jac__)
            for grp in fitgroups.__grps__:            
                if active_only and not grp.__active__: continue
                STATS = grp.calculate_outlier_statistics()
            for stat in STATS:
                fitgroups.collect(stat,active_only=active_only)
        
        outlier_stats.append(['leverage',fitgroups.__leverage__])
        outlier_stats.append(['student',fitgroups.__student__])
        outlier_stats.append(['dffits',fitgroups.__DFFITS__])
        outlier_stats.append(['cook',fitgroups.__cook__])

        
    # prepare header
    header = list(input_names) + [output_column] + ['obs','calc','weight','wresid','uresid'] + \
        [stat_name for stat_name,_ in outlier_stats]
        
    col = j.Collection(); col.order = header
    
    #print('header>>>',header)
        
    # loop through vals to create a full stat collection
    for i in range(length):
        
        input_vals = input_matrix[i]
        
        item = {
            'obs': output[i],
            'calc': model_calc[i],
            'weight': full_weights[i],
            'wresid': weighted_residuals[i],
            'uresid': unweighted_residuals[i],
        }
        
        item[output_column] = output[i]
        
        for input_name,input_val in zip(input_names,input_vals):
            item[input_name] = input_val
        
        for stat_name,stat_vals in outlier_stats:
            item[stat_name] = stat_vals[i]
            
        col.update(item)
        
    #print(col.keys())
        
    stat_buffer = col.tabulate(raw=True)
    
    # tabulate collection to output, if specified
    if stream:
        stream.write('=======================================\n')
        stream.write('RESIDUAL STATISTICS\n')
        stream.write('=======================================\n')
        stream.write(stat_buffer)
        
    if save:
        # save outlier statistics to file
        col.export_csv('%s.resids.csv'%model_name)
        
        # save Pareto front based on outlier statistics (if provided)
        if outlier_stats:
        
            print('Calculating Pareto set for the Leverage/Cook pairs...')
        
            def func_factory():
                count = -1
                def get_id(v):
                    nonlocal count
                    count += 1
                    return count
                return get_id
        
            leverage,cook,ids = col.getcols(['leverage','cook','id'],
                functions={'id':func_factory()})
        
            costs = np.array((leverage,cook))
            costs = costs.T
            mask = is_pareto_efficient(costs,return_mask=False)
            ids_pareto = list(np.array(ids)[mask])
            col_pareto = col.subset(ids_pareto)
            col_pareto.order = col.order
            pareto_file = '%s.pareto.csv'%model_name
            col_pareto.export_csv(pareto_file)
            
            print('...saved to %s'%pareto_file)

def calculate_and_save_statistics(CONFIG,model,fitgroups,stream=sys.stdout,active_only=True):
    fitgroups.calculate(model,active_only=active_only)
    #fitgroups.collect_calc_vals(active_only=active_only)
    #fitgroups.collect_resids(active_only=active_only)
    calculate_fitpar_statistics(CONFIG,fitgroups,save=True,stream=stream,active_only=active_only)
    calculate_group_statistics(CONFIG,fitgroups,save=True,stream=stream,active_only=active_only)
    calculate_residual_statistics(CONFIG,fitgroups,save=True,stream=stream,active_only=active_only)

def stat(CONFIG):
    """ Calculate fit statistics """
    
    # Get options from the config file.
    model_name = CONFIG['MODEL']['model']
    stat_file = CONFIG['STAT']['stat_file']
    output_symbolic_func = to_bool(CONFIG['STAT']['output_symbolic_func'])
    
    # Get model and read parameters from file.
    model = load_model(CONFIG)
    
    # Load fitgroups
    fitgroups,_ = read_fitgroups(CONFIG,verbose=1)
    
    # Load Jacobian
    jacfile = model.__class__.__name__+'.jac.npy'
    if os.path.isfile(jacfile):
        jac = np.load(jacfile)
    else:
        jac = None
        
    fitgroups.__jac__ = jac
                       
    # Print symbolic function.
    if output_symbolic_func:
        print('\n\nSYMBOLIC FUNC (ORIGINAL)=====>')
        print(f.model_final.__symbolic_func__)
        print('\n\nSYMBOLIC FUNC (SUBS)=====>')
        dct = {psym.get_value():f.model_final.__params__[psym.__name__].get_value() \
            for psym in f.model_final.__symbolic_params__} # parameters substitution
        expr_sub = f.model_final.__symbolic_func__.subs(dct) # substitute parameters and inputs 
        print(expr_sub)
    
    # Save statistics.
    with open(stat_file,'w') as stat_stream:
        calculate_and_save_statistics(CONFIG,model,fitgroups,stream=stat_stream,active_only=True)

#################
#### CODEGEN ####
#################

def get_arguments(CONFIG):
    """ Get argument names of the model """    
    
    arguments_line = CONFIG['MODEL']['arguments']
    
    if not arguments_line:
        print('ERROR: MODEL.arguments list should not be empty!')
        sys.exit()
        
    argnames = [argname.strip() for argname in arguments_line.split(';')]
    if not argnames[-1]: argnames = argnames[:-1]
    
    return argnames
    
def parse_gridspec(CONFIG,gridspec_line):
    """ Parse the grid specifications for plotting and compariong with generated code. 
        Format of the gridspec line: X=XMIN:XSTEP:XMAX; Y=YMIN:YSTEP:YMAX; ... """
        
    if not gridspec_line:
        print('ERROR: gridspec must be set before using plotting and codegen.')
        sys.exit()
        
    arguments = get_arguments(CONFIG)
    
    argpos = {} # indexes of arguments for assign the indexes_unfixed var
    for i,argname in enumerate(arguments):
        argpos[argname] = i
    
    gridspec_dict = {}
    indexes_unfixed = []
    argset = set()
    for grsp in gridspec_line.split(';'):
        grsp = grsp.strip()
        argname,grstr = grsp.split('=')
        grstr_split = grstr.split(':')
        if len(grstr_split)==3:
            cmin,cstep,cmax = grstr_split
            cmin = float(cmin)
            cstep = float(cstep)
            cmax = float(cmax)
            if cmin+cstep>cmax: cstep = cmax-cmin
            if cstep==0:
                value = np.array([cmin])
            else:
                value = np.arange(cmin,cmax,cstep)
            indexes_unfixed.append(argpos[argname])
        elif len(grstr_split)==2:
            cmin,cmax = grstr_split
            cmin = float(cmin)
            cmax = float(cmax)
            if cmin==cmax:
                value = np.array([cmin])
            else:
                value = np.linspace(cmin,cmax,100)
            indexes_unfixed.append(argpos[argname])
        elif len(grstr_split)==1:
            value = float(grstr_split[0])
        else:
            print('ERROR: bad gridspec format: %s'%gridspec_line)
            sys.exit()
        gridspec_dict[argname] = value
        argset.add(argname)
    
    remainder = set(arguments)-argset
    if remainder:
        print('ERROR: %s args not in gridspec'%remainder)
        sys.exit()
        
    gridspec = []
    for i,argname in enumerate(arguments):
        val = gridspec_dict[argname]
        if type(val) is np.ndarray:
            gridspec.append((argname,val))
        elif type(val) is float:
            gridspec.append((argname,[val]))
        else:
            raise Exception('unsupported type: %s'%type(val))
    
    return gridspec,indexes_unfixed

def codegen(CONFIG):
    """ Code generator for the fit model. """
    
    # Get options from the config file.    
    create_fortran = to_bool( CONFIG['CODEGEN']['create_fortran'] )
    compare_fortran = to_bool( CONFIG['CODEGEN']['compare_fortran'] )
    compiler = CONFIG['CODEGEN']['compiler_fortran']
    
    # =====================
    # GENERATE FORTRAN CODE
    # =====================

    # Get model and read parameters from file.
    model = load_model(CONFIG)

    # Get grid to calculate the generated routines on.
    gridspec,_ = parse_gridspec(CONFIG,CONFIG['CODEGEN']['gridspec'])
    grid_compfort = Grid(*reduce(lambda x,y:x+y,gridspec))
    calc_model_python = model.calculate(grid_compfort)
        
    # Create compiled Fortran gridcalc
    if create_fortran:
        model.generate_fortran_90()
        create_pes_interface(model)
        create_pots_interface()

    # Compare Python version with fortran is needed.
    if compare_fortran:
        model.compile_fortran_90(compiler)
        calc_model_fortran = model.__function_fortran__(grid_compfort)
        calc_model_fortran = np.squeeze(calc_model_fortran) # fixes bug when some dimensions are "degenerate"
        abs_err_compfort = np.abs(calc_model_python - calc_model_fortran)
        abs_err_compfort_min = np.min(abs_err_compfort)
        abs_err_compfort_max = np.max(abs_err_compfort)
        COMPFORT_THRESH = 1.0E-7
        print('\n\n======================================================================')
        print('VALUES FOR PYTHON-FORTRAN COMPARISON: MIN=%E, MAX=%E, THRESH=%E'%\
            (abs_err_compfort_min,abs_err_compfort_max,COMPFORT_THRESH))
        if abs_err_compfort_max>COMPFORT_THRESH:
            print('WARNING: maximal discrepancy between Python and generated Fortran models'
                'is too large: %E'%abs_err_compfort_max)
        else:
            print(' ...COMPARISON TEST PASSED!')
        print('======================================================================')
        print('\n')

##############
#### PLOT ####
##############

def plot_residuals(CONFIG):
    """ Plot fit residuals versus a given coordinate. """
    
    # Get options from the config file.
    resids_weighted = to_bool(CONFIG['PLOTTING']['resids_weighted'])
    #module_name = CONFIG['MODEL']['model']
    #model_name = module_name
    #modfile = model_name+'.model'

    # Import model module (in other case, deserialization is not working).
    #module_path = os.path.join('./',module_name+'.py') # relative path
    #module = import_module_by_path(module_name,module_path)

    # Get input and output names
    INPUTS = parse_input_columns(CONFIG)
    OUTPUT = CONFIG['DATA']['output_column']
    dim = len(INPUTS)
    
    # Get addiitonal model list to compare with
    compare_with_models = CONFIG['PLOTTING']['compare_with_models'].strip()

    # X axis for residuals.
    resids_x_axes = CONFIG['PLOTTING']['resids_x_axes']   
    resids_x_axes = resids_x_axes.strip()
    if not resids_x_axes: resids_x_axes=OUTPUT

    # Get model and read parameters from file.
    #model = deserialize_model(modfile)
    model = load_model(CONFIG)
    
    # Get data collection.
    _,col_data = read_fitgroups(CONFIG,verbose=1)
        
    # Prepare observable data.
    ids_sort = col_data.sort(resids_x_axes)
    DATA = col_data.getcols(
        #INPUTS+[OUTPUT,resids_x_axes,'N'],
        INPUTS+[OUTPUT,resids_x_axes],
        IDs=ids_sort)
    input_columns = DATA[:dim]
    output_column, = DATA[dim:dim+1]
    #resid_axis_data,nn = DATA[dim+1:]
    resid_axis_data = DATA[dim+1]
    grid = List(*reduce(lambda x,y:x+y,zip(INPUTS,input_columns))) # list-based grid
    
    # Do the plotting part.
    #model_vals = model.calculate(grid)
    
    models = [model]
    if compare_with_models:
        model_names = parse_semicolon_list(CONFIG,'PLOTTING','compare_with_models')
        for model_name in model_names:
            models.append(load_model_rec(model_name))
    
    # =============
    # UPPER SUBPLOT
    # =============
    
    ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3)
    ax1.get_xaxis().get_major_formatter().set_useOffset(False) # set normal format

    leg = ['DATA']
    
    # plot observable data
    #plt.plot(resid_axis_data,output_column,'ro'); leg.append('data')
    plt.scatter(resid_axis_data, output_column, fc='none', ec='r', s=100)

    model_calc_vals = []
    for model,color,marker in zip(models,cycle(COLORS),cycle(MARKERS)):

        model_vals = model.calculate(grid)
        plt.plot(resid_axis_data,model_vals,marker,color=color)
        leg.append(model.__class__.__name__)
        model_calc_vals.append(model_vals)

        #if SHOW_N:
        #    for n,v,raxis in zip(nn,output_column,resid_axis_data):
        #        ax1.text(raxis,v,'%d'%n)

    plt.title('residuals')
    plt.legend(leg)
    plt.grid(True)
    plt.ylabel(OUTPUT)
    
    # =============
    # LOWER SUBPLOT
    # =============

    plt.subplot2grid((4,1),(3,0),rowspan=1,sharex=ax1)

    #leg = []
    for model_vals,model,color,marker in zip(model_calc_vals,models,cycle(COLORS),cycle(MARKERS)):
        plt.plot(resid_axis_data,output_column-model_vals,marker+'-',color=color)
        #leg.append(model.__class__.__name__)
    
    plt.xlabel(resids_x_axes)
    plt.ylabel('OBS-CALC')
    #plt.legend(leg)
    plt.grid(True)
    
    plt.show()    

def get_argument_bindings(CONFIG):
    """
    Parse coordinate bindings to the names of the data columns.
    Bindings should have the following format:
    coord1:datacol1;coord2:datacol2;...
    Data column names should correspond to the ones mentioned 
    in the DATA section.
    Order of the binding MUST correspond to the 
    argument of the model's __func__ method. 
    """
    model_arguments = CONFIG['MODEL']['arguments']
    input_columns = CONFIG['DATA']['input_columns']
    bindings = []
    for pair in zip(model_arguments.split(';'),input_columns.split(';')):
        if not pair: continue
        argname = pair[0].strip()
        colname = pair[1].strip()
        bindings.append((argname,colname))
    return bindings

def make_title(gridspec,indexes_unfixed,bindings):
    """ Make title for all plotting functions """
    gridspec_dict = dict(gridspec)    
    n_args = len(bindings)
    indexes_fixed = list(set(range(n_args))-set(indexes_unfixed))
    pieces = []
    for index in indexes_fixed:
        argname = bindings[index][0]
        val = gridspec_dict[argname][0]
        pieces.append('%s=%.3f'%(argname,val))
    title = '; '.join(pieces)
    return title

def get_model_components(CONFIG):
    """ Get the list of the model components separated by semicolon."""
    components_line = CONFIG['PLOTTING']['model_components']
    components = []
    for compname in components_line.split(';'):
        compname = compname.strip()
        if not compname: continue
        components.append(compname)
    return components

def plot_sections(CONFIG):
    """ Plot model cuts and compare to the data points (if any). """
    
    # Get argument names.
    argnames = get_arguments(CONFIG)
    n_args = len(argnames)
    
    # Get coordinate bindings.
    bindings = get_argument_bindings(CONFIG)
    bindings_dict = dict(bindings)
    
    # Get main model name
    main_model_name = CONFIG['MODEL']['model']
    
    # Get gridspec and indexes of unfixed arguments
    gridspec,indexes_unfixed = \
        parse_gridspec(CONFIG,CONFIG['PLOTTING']['gridspec'])
        
    n_unfixed = len(indexes_unfixed)
        
    gridspec_dict = dict(gridspec)
        
    indexes_fixed = list(set(range(n_args))-set(indexes_unfixed))
    
    scatter_opacity = to_float( CONFIG['PLOTTING']['scatter_opacity'] )
    marker_size = to_float( CONFIG['PLOTTING']['marker_size'] )
    #resize_by_weights = to_bool( CONFIG['PLOTTING']['resize_by_weights'] )
    surface_opacity = to_float( CONFIG['PLOTTING']['surface_opacity'] )
    calculate_components = to_bool( CONFIG['PLOTTING']['calculate_components'] )
    outlier_stats_type = CONFIG['STAT']['outlier_stats_type']
    plot_outlier_stats = to_bool( CONFIG['PLOTTING']['plot_outlier_stats'] )
    exclude = CONFIG['DATA']['exclude']
    annotate_field = CONFIG['PLOTTING']['annotate_field']

    # Get input and output names
    INPUTS = parse_input_columns(CONFIG)
    OUTPUT = CONFIG['DATA']['output_column']
    
    # Get addiitonal model list to compare with
    compare_with_models = CONFIG['PLOTTING']['compare_with_models'].strip()
    
    # Import model module (in other case, deserialization is not working).   
    model = load_model(CONFIG,verbose=0)
    
    # Get fitgroups and calculate statistics.
    fitgroups,col_data = read_fitgroups(CONFIG,verbose=1)
    
    # Get model components list.
    components = get_model_components(CONFIG)
    
    # Do plotting stuff.
    if n_unfixed==2:
        fig = plt.figure()
        plotter = fig.add_subplot(111, projection='3d')
    elif n_unfixed==1:
        plotter = plt
    else:
        raise NotImplementedError

    print('')

    leg = []
    
    # Calculate outlier statistics if needed.
    if plot_outlier_stats:
        jacfile = model.__class__.__name__+'.jac.npy'
        print('Loading Jacobian matrix from %s.'%jacfile)
        fitgroups.__jac__ = np.load(jacfile)
        print('Calculating outlier statistics...')
        fitgroups.calculate(model,active_only=True)
        with open(os.devnull,'w') as devnull:
            calculate_residual_statistics(CONFIG,fitgroups,save=False,
                stream=devnull,active_only=True)
        print('...DONE')

    # Plot data points.
    for fitgroup,grp_color in zip(fitgroups,cycle(COLORS)):
        
        # find index of all data points whose fixed columns correspond to gridspec
        ind = np.full(fitgroup.__inputgrid__.__shape__,True)
        for index_fixed in indexes_fixed:
            argname_fixed,colname_fixed = bindings[index_fixed]
            val_fixed = gridspec_dict[argname_fixed][0]
            ind *= abs(fitgroup[colname_fixed]-val_fixed)<0.001
            
        n_tot = len(ind)
        n_sel = np.count_nonzero(ind)
        
        if n_sel==0: continue
        
        leg.append('%s (%d/%d)'%(fitgroup.__name__,n_sel,n_tot))
        
        meshes = fitgroup.__inputgrid__.get_meshes(ind)
        plot_data = [meshes[i] for i in indexes_unfixed]
        output_data = fitgroup.__output__[ind]
        if plot_outlier_stats:
            stat = fitgroup.stats(outlier_stats_type)[ind]
            stat_title = '(%s)'%outlier_stats_type
            stat_colors = stat
        else:
            stat_title = ''
            stat_colors = grp_color
        
        if n_unfixed in [1,2]:
            sc = plotter.scatter(*plot_data,output_data,c=stat_colors,
                        s=marker_size,alpha=scatter_opacity)
        else:
            raise NotImplementedError

    # Define condition for finding 1D section.
    def section_cond(v):
        cond = True
        for index_fixed in indexes_fixed:
            argname_fixed,colname_fixed = bindings[index_fixed]
            val_fixed = gridspec_dict[argname_fixed][0]
            cond_ = v[colname_fixed]==val_fixed
            cond = cond and cond_
        return cond
            
    # Plot excluded data points.
    if exclude:
        col_exclude = j.import_csv(exclude)
        col_exclude = col_exclude.subset(col_exclude.ids(lambda v: OUTPUT in v))
        col_exclude = col_exclude.subset(col_exclude.ids(section_cond))
        plot_data = col_exclude.getcols(
            [bindings[index_unfixed][1] \
                for index_unfixed in indexes_unfixed]+[OUTPUT]
        )
        if plot_data[0]:
            plotter.scatter(*plot_data,s=marker_size+50,
                facecolors='none', color='red',
                alpha=scatter_opacity)
            leg.append('data excluded')
            
    # Plot annotations.
    if annotate_field:
        col_data = col_data.subset(col_data.ids(section_cond))
        plot_data = col_data.getcols(
            [bindings[index_unfixed][1] \
                for index_unfixed in indexes_unfixed]+[OUTPUT]+[annotate_field]
        )
        for d in zip(*plot_data):
            plotter.text(*d)
    
    # Plot models.
    g = Grid(*reduce(lambda x,y:x+y,gridspec))
    meshes = g.get_meshes()
    calc_data = [meshes[i] for i in indexes_unfixed]
    
    models = [model]
    model_names = [main_model_name]
    if compare_with_models:
        aux_model_names = parse_semicolon_list(CONFIG,'PLOTTING','compare_with_models')
        model_names += aux_model_names
        for model_name in aux_model_names:
            models.append(load_model_rec(model_name))
        
    for model,color,model_name in zip(models,cycle(COLORS),model_names):
            
        if calculate_components:
            results = model.calculate_components(g,compnames=components)
        else:
            results = model.calculate(g); results = [results]
            
        for res,compname in zip(results,['%s (FULL MODEL)'%model_name]+\
                ['%s (%s)'%(model_name,c) for c in components]):
            print('\nplotting %s (%s) for %s'%(compname,color,model_name))
            leg.append(compname)
            if n_unfixed==2:
                surf = plotter.plot_surface(*calc_data,res,alpha=surface_opacity,color=color)
                # solution for Matplotlib >=3.3
                # https://stackoverflow.com/questions/4536103/how-can-i-upgrade-specific-packages-using-pip-and-a-requirements-file
                surf._facecolors2d = surf._facecolor3d
                surf._edgecolors2d = surf._edgecolor3d
            elif n_unfixed==1:
                plotter.plot(*calc_data,res,color=color)
            else:
                raise NotImplementedError

    plt.title(make_title(gridspec,indexes_unfixed,bindings))
    
    if n_unfixed==2:
        plotter.set_xlabel('%s'%bindings_dict[argnames[indexes_unfixed[0]]])
        plotter.set_ylabel('%s'%bindings_dict[argnames[indexes_unfixed[1]]])
        plotter.set_zlabel('%s'%OUTPUT)
    elif n_unfixed==1:
        plotter.xlabel('%s'%bindings_dict[argnames[indexes_unfixed[0]]])
        plotter.ylabel('%s'%OUTPUT)
    else:
        raise NotImplementedError
    
    if plot_outlier_stats: 
        plt.colorbar(sc)

    #leg += ['calc_model']
    plotter.legend(leg)

    plt.grid(True)
    plt.show()

def plot_fit_history(CONFIG):
    """ Plot fitting history for each data group. """

    # Get options from the config file.
    model_name = CONFIG['MODEL']['model']
    fit_history_stat = CONFIG['PLOTTING']['fit_history_stat']
    fit_history_logscale = to_bool( CONFIG['PLOTTING']['fit_history_logscale'] )
    
    # Open history file in JSON format.
    with open(model_name+'.fit_history') as f:
        fit_history = json.load(f)
            
    # Start plotting.
    GROUPS = {}
    
    for niter,global_item in enumerate(fit_history):
        stat_dicthash = global_item['stat_dicthash']
        for key in stat_dicthash:
            group_item = stat_dicthash[key]
            group_name = group_item['GROUP']
            val = group_item[fit_history_stat]
            if group_name not in GROUPS:
                GROUPS[group_name] = [(niter,val)]
            else:
                GROUPS[group_name].append((niter,val))
    
    leg = []
    
    for group_name,marker,color in zip(
        list(GROUPS.keys()),
        cycle(MARKERS),
        cycle(COLORS)
        ):
        
        pairs = GROUPS[group_name]
        #pairs = sorted(pairs)
        niters,vals = list(zip(*pairs))
        pl.plot(niters,vals,marker+'-',color=color)
        leg.append(group_name)
    
    pl.xlabel('N_iter')
    pl.ylabel(fit_history_stat)
    pl.title('History for '+fit_history_stat)
    pl.legend(leg)
    pl.grid(True)
    if fit_history_logscale: pl.yscale('log')
    pl.show()

def plot_pareto(CONFIG):
    """ Plot Pareto set for the most influential outliers. """
    
    # Get options from the config file.
    model_name = CONFIG['MODEL']['model']

    # Open cached file with residuals Pareto set and plot.
    col = j.import_csv('%s.resids.csv'%model_name)
    col_pareto = j.import_csv('%s.pareto.csv'%model_name)
    
    pl.plot(*col.getcols(['leverage','cook']),'bo')
    pl.plot(*col_pareto.getcols(['leverage','cook']),'ro')
    pl.xscale('log')
    pl.yscale('log')
    pl.xlabel('leverage')
    pl.ylabel('cook')
    pl.show()

def plot(CONFIG):
    """ interface for the plotting section. Supports a number of standard plots. """
    
    # Get options from the config file.
    plot_mode = CONFIG['PLOTTING']['plot_mode']
    
    # Switch on the plotting options.
    if plot_mode=='residuals':
        plot_residuals(CONFIG)
    elif plot_mode=='sections':
        plot_sections(CONFIG)
    elif plot_mode=='history':
        plot_fit_history(CONFIG)
    elif plot_mode=='pareto':
        plot_pareto(CONFIG)
    else:
        print('ERROR: unknown plotting mode "%s"'%plot_mode)
        sys.exit()

def plot_multicut(CONFIG):
    """ plot multiple 1D cuts and multiple models """

    # Get coordinate bindings.
    bindings = get_argument_bindings(CONFIG)
    bindings_dict = dict(bindings)
    
    # Get options from config.
    exclude = CONFIG['DATA']['exclude']
    E_COLUMN = CONFIG['DATA']['output_column']
    show_legend = to_bool(CONFIG['MULTICUT']['show_legend'])
    unfixed_arg_name = CONFIG['MULTICUT']['argument']
    plot_model = to_bool(CONFIG['MULTICUT']['plot_model'])
    fullgrid_calc = to_bool(CONFIG['MULTICUT']['fullgrid_calc'])
    show_lines = to_bool(CONFIG['MULTICUT']['show_lines'])
    
    all_arg_names = [c[0] for c in bindings]
    fixed_arg_names = [c[0] for c in bindings if c[0]!=unfixed_arg_name]

    # Annotations.
    annotate_field = CONFIG['MULTICUT']['annotate_field'].strip()
    if not annotate_field: 
        plot_annotations = False
        annotate_field = '__ID__' # stub
    else:
        plot_annotations = True
        
    # Input names.
    all_input_names = [c[1] for c in bindings]
    unfixed_input_name = bindings_dict[unfixed_arg_name]
    fixed_input_names = [bindings_dict[c] for c in fixed_arg_names]

    # Parse gridspec. 
    gridspec,indexes_unfixed = \
        parse_gridspec(CONFIG,CONFIG['MULTICUT']['gridspec'])
    gridspec_dict = dict(gridspec)
        
    # Get array-like options.    
    bnds_dict = {
        bindings_dict[c]: [gridspec_dict[c][0],gridspec_dict[c][-1]] \
            for c in bindings_dict
    }
    
    # Get main model name
    main_model_name = CONFIG['MODEL']['model']

    # Get addiitonal model list to compare with
    compare_with_models = CONFIG['MULTICUT']['compare_with_models']

    # Import model module (in other case, deserialization is not working).   
    main_model = load_model(CONFIG)
    
    # Get fitgroups and calculate statistics.
    fitgroups,col_data = read_fitgroups(CONFIG,verbose=1,exclude=False)
    
    # Get excluded points separately
    if exclude: col_excluded = j.import_csv(exclude)

    # Filter data points to account for bounds on coordinates.
    col_data = col_data.subset(col_data.ids(
        lambda v: reduce(lambda a,b: a and b,[
            bnds_dict[c][0] <= v[c] <= bnds_dict[c][1] \
                for c in bnds_dict
        ])
    ))

    # Define legend line.
    def get_legend_line(fixed_input_names,fixed_input_vals,postfix):
        return ', '.join(
            ['%s=%.3f'%(name,val) for name,val in zip(
                fixed_input_names,fixed_input_vals)]) + ' %s'%postfix

    # Start plotting.
    fig = plt.figure()
    ax = fig.add_subplot()

    excluded_points = []
    if exclude:
        col_excluded_ID_vals = col_excluded.group(
            lambda v:tuple(v[c] for c in all_input_names))
    else:
        col_excluded_ID_vals = set()

    leg = []
    color_gen = cycle(COLORS)
    marker = 'o'
    data_linestyle = '-' if show_lines else ''
            
    grpi_sections = col_data.group(lambda v: tuple(v[c] for c in fixed_input_names))
    grpi_sections_keys = sorted(grpi_sections)
            
    for fixed_input_vals in grpi_sections_keys:
                    
            col = col_data.subset(grpi_sections[fixed_input_vals])
                    
            annot_vals,unfixed_input_vals,E_vals = col.getcols(
                [annotate_field,unfixed_input_name,E_COLUMN],
                IDs=col.sort(unfixed_input_name),
            )

            color = next(color_gen)

            ax.plot(unfixed_input_vals,E_vals,marker=marker,color=color,
                linestyle=data_linestyle)
            
            lookup = {coord_name:coord_val \
                for coord_name,coord_val in zip(fixed_input_names,fixed_input_vals)}
            
            for xx,ee,aa in zip(unfixed_input_vals,E_vals,annot_vals):
                # plot grid points IDs
                if plot_annotations: 
                    ax.text(xx,ee,aa)
                # collect exclude points
                lookup[unfixed_input_name] = xx
                key = tuple(lookup[c] for c in all_input_names)
                if key in col_excluded_ID_vals:
                    excluded_points.append((xx,ee))
            
            leg.append(get_legend_line(fixed_input_names,fixed_input_vals,'(points)'))
        
    indexes_unfixed_ = [all_arg_names.index(unfixed_arg_name)] # specially for plotting models
        
    def plot_models_on_grid(g,fixed_input_names,fixed_input_vals,color,marker_flag=False,append_legend=True):
            
        meshes = g.get_meshes()
        calc_data = [meshes[i] for i in indexes_unfixed_]
    
        models = [main_model]
        model_names = [main_model_name]
        if compare_with_models:
            aux_model_names = parse_semicolon_list(CONFIG,'MULTICUT','compare_with_models')
            model_names += aux_model_names
            for model_name in aux_model_names:
                models.append(load_model_rec(model_name))
        
        for model,marker,model_name in zip(models,cycle(MARKERS),model_names):
            results = model.calculate(g)            
            legline = get_legend_line(fixed_input_names,fixed_input_vals,'(%s)'%model_name)
            print('calculating',legline)
            if append_legend:
                leg.append(legline)
            if not marker_flag:
                marker = None
            ax.plot(*calc_data,results,color=color,marker=marker,linestyle='--')
        
    # Plot models.
    
    #marker_flag = fullgrid_calc
    marker_flag = False
    
    if plot_model:
        
        # If fullgrid_calc=True, plot model on full grid.        
        if fullgrid_calc:
            
            print('PLOTTING FULL GRID MODEL CUTS')
            
            color_gen = cycle(COLORS)
            
            fixed_arg_vals = [gridspec_dict[arg] for arg in fixed_arg_names]
            fixed_arg_vals_meshgrid = np.meshgrid(*fixed_arg_vals)
            fixed_arg_vals_meshgrid = [c.flatten() for c in fixed_arg_vals_meshgrid]
            col_fixed_fullgrid = j.Collection()
            col_fixed_fullgrid.update([
                {arg:val for arg,val in zip(fixed_arg_names,tup) } \
                    for tup in zip(*fixed_arg_vals_meshgrid)
            ])
            
            print('Fixed values for the full grid calculation:')
            col_fixed_fullgrid.tabulate(fixed_arg_names)
            
            grpi_sections_ = col_fixed_fullgrid.group(
                lambda v: tuple(v[c] for c in fixed_arg_names))
            grpi_sections_keys_ = sorted(grpi_sections_)
            
            for fixed_arg_vals in grpi_sections_keys_:
                
                lookup = {coord_name:coord_val \
                    for coord_name,coord_val in zip(fixed_arg_names,fixed_arg_vals)}
                    
                gridspec_ = list(
                    (c,gridspec_dict[c]) \
                        if c==unfixed_arg_name else (c,[lookup[c]]) \
                    for c in all_arg_names)
                    
                g = Grid(*reduce(lambda x,y:x+y,gridspec_))
            
                color = next(color_gen)
                            
                plot_models_on_grid(g,fixed_arg_names,fixed_arg_vals,color,
                    marker_flag,append_legend=False)
        
        color_gen = cycle(COLORS)
        
        # Plot model on sections where there are data.
        print('PLOTTING MODEL VS DATA CUTS')
        for fixed_input_vals in grpi_sections_keys:
            
            lookup = {coord_name:coord_val \
                for coord_name,coord_val in zip(fixed_arg_names,fixed_input_vals)}
            
            gridspec_ = list(
                (c,gridspec_dict[c]) \
                    if c==unfixed_arg_name else (c,[lookup[c]]) \
                for c in all_arg_names)
                
            g = Grid(*reduce(lambda x,y:x+y,gridspec_))
            
            color = next(color_gen)
            
            plot_models_on_grid(g,fixed_arg_names,fixed_input_vals,color,
                marker_flag,append_legend=not fullgrid_calc)
                                
    if excluded_points:
        # plot excluded points
        ax.scatter(*zip(*excluded_points), s=150, facecolors='none', color='black')
        leg.append('bad points')

    ax.set_xlabel(unfixed_arg_name)
    ax.set_ylabel(E_COLUMN)

    #plt.title(filename)

    if show_legend: 
        plt.legend(leg)

    plt.show()
    
def calc(CONFIG):
    """ Calculate fitted model on a grid. """
    
    # Get options from the config file.    
    output_file = CONFIG['CALC']['output_file']
    
    # Get input and output names.
    inputs = parse_input_columns(CONFIG)
    output = CONFIG['DATA']['output_column']

    # Load model.
    model = load_model(CONFIG)
    
    # Get gridspec and grid.
    gridspec,_ = \
        parse_gridspec(CONFIG,CONFIG['CALC']['gridspec'])
    grid_calc = Grid(*reduce(lambda x,y:x+y,gridspec))

    # Calculate model and flatten meshes.
    calc_model = model.calculate(grid_calc)
    calc_model = calc_model.flatten()
    meshes = grid_calc.get_meshes(flat=True)
        
    # Create and write output CSV file.
    with open(output_file,'w') as f:
        # header
        f.write(';'.join(inputs+[output])+'\n')
        # body
        for vals in zip(*(meshes+[calc_model])):
            buf = ';'.join([str(val) for val in vals])+'\n'
            f.write(buf)
            
    print('Calculated model values are saved to %s'%output_file)
