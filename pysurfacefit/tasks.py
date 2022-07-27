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

# type conversions minding the None value
to_type = lambda val,typ: typ(val) if val not in (None,'') else None
to_int = lambda val: to_type(val,int)
to_float = lambda val: to_type(val,float)
to_bool = lambda val: val.lower() in ("yes", "true", "t", "1")

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
    with open(os.path.join(project_dir,filename),'w') as f:
        f.write(content)
    return filename

def startproject(CONFIG):
    project_dir = CONFIG['GENERAL']['project']
    create_dir(project_dir)
    config_name = 'config.ini'
    config_path = os.path.join(project_dir,config_name)
    CONFIG.save(config_path)
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

def parse_dataspec(CONFIG):
    filename = CONFIG['DATA']['dataspec']
    col_dataspec = j.import_fixcol(filename)
    return col_dataspec

def read_fitgroups(CONFIG,verbose=False):
    """ Read data fitgroups based on the config file information. """
    
    # Create unified collection for plotting.
    col_data = j.Collection() # unified collection for plotting

    # Get the config options.
    datafile = CONFIG['DATA']['datafile']

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
    
    # First, accumulate all datafiles in one place
    col_dataspec = j.Collection()
    if datafile:
        col_dataspec.update({
            'alias':'default','path':datafile,
            'wht_mul':1.0,'type':None,'include':1})
    
    # Second, get fitgroups from the data in the csv file.
    col = parse_dataspec(CONFIG)
    col_dataspec.update(col.getitems())
    
    GROUPS_DATA = []
    for item in col_dataspec.getitems():
        # Get fitgroup properties
        fitgroup_path = item['path'].strip()
        fitgroup_name = item['alias'].strip()
        fitgroup_type = item['type']
        fitgroup_wht_mul = item['wht_mul']
        if verbose: print('\n===============================')
        if verbose: print(fitgroup_name)
        if verbose: print('===============================')
        # Get points from the CSV file
        col = j.Collection(); col.import_csv(fitgroup_path)
        # Apply data cutoff
        col = col.subset(col.ids(lambda v: GLOBAL_CUTOFF_MIN<=v[OUTPUT]<=GLOBAL_CUTOFF_MAX))
        # Get the data columns 
        input_columns = col.getcols(INPUTS)
        output_column = col.getcol(OUTPUT)
        # Calculate data weights 
        col.assign('fit_weight',wht_fun)
        data_whts = col.getcol('fit_weight')
        # Print collection contents
        if verbose:
            col.tabulate(['N',*(INPUTS+[OUTPUT]),'fit_weight'])
        # Prepare fitting grid (list)
        grid = List(*reduce(lambda x,y:x+y,zip(INPUTS,input_columns))) # list-based grid
        # Prepare fitgroup for points
        if fitgroup_type is None:
            FitPointsType = FitPoints
        elif fitgroup_type==-1:
            FitPointsType = PenaltyPointsDown
        elif fitgroup_type==1:
            FitPointsType = PenaltyPointsUp
        else:
            print('ERROR: unknown dataspec type: %s'%str(fitgroup_type))
            sys.exit()
        print('Treating %s fitgroup as %s'%(fitgroup_name,FitPointsType.__name__))
        GRP_DATA = FitPoints(
            name=fitgroup_name,input=grid,output=output_column,
            whtmul=fitgroup_wht_mul,whts=data_whts,active=True
        )
        # Add to the plotting collection
        col_data.update(col.getitems())
        # Append to groups of data
        GROUPS_DATA.append(GRP_DATA)
        
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

def load_model(CONFIG):
    """ Load model from its module.
        Model should read its parameters automatically! """
    module_name = CONFIG['MODEL']['model']
    #module_path = os.path.abspath(module_name+'.py') # absolute path
    module_path = os.path.join('./',module_name+'.py') # relative path
    module = import_module_by_path(module_name,module_path)
    model_object = 'model'
    model = getattr(module,model_object)
    return model

def fit(CONFIG):
    """ Start fitting procedure. """
    
    # Get options from the config file.
    weighted_fit = CONFIG['FIT']['weighted_fit']
    rubber_on = CONFIG['FIT']['rubber_on']
    fitting_method = CONFIG['FIT']['fitting_method']
    analytic_jacobian = CONFIG['FIT']['analytic_jacobian']
    #module_name = CONFIG['MODEL']['model']
    #module_path = os.path.abspath(module_name+'.py') # absolute path
    #module_path = os.path.join('./',module_name+'.py') # relative path
    fitopts = parse_fit_options(CONFIG)
    #model_name = 'Model'
    #model_name = module_name
    #model_object = 'model'
    #model = getattr(__import__(module_name),model_name)
    #module = import_module_by_path(module_name,module_path)
    #model = getattr(module,model_object)
    model = load_model(CONFIG)
    fitgroups,_ = read_fitgroups(CONFIG)
    #fitfile = module_name+'.fit'
    #modelfile = module_name+'.model'
    
    # Create Fitter object from scratch.
    f = Fitter(model=model,fitgroups=fitgroups,
            weighted=weighted_fit,rubber_on=rubber_on,
            method=fitting_method,jac=analytic_jacobian,
            **fitopts);
    
    # If fit file is absent, create a fit from scratch.
    #if not os.path.exists(fitfile):
    #    f = Fitter(model=model,fitgroups=fitgroups,
    #        weighted=weighted_fit,rubber_on=rubber_on,
    #        method=fitting_method,jac=analytic_jacobian,
    #        **fitopts);
    #else:
    #    f = deserialize_fit(fitfile)
    #    f.__fitgroups__ = fitgroups    
        
    # Start fit.
    f.fit()
    
    # Save fit.
    #serialize_fit(f,fitfile)

    # Save model.
    #serialize_model(f.model_final,modelfile)
    f.model_final.save_params(model.__class__.__name__+'.csv')
    
##############
#### STAT ####
##############
    
def stat(CONFIG):
    """ Calculate fit statustics """
    
    # Get options from the config file.
    module_name = CONFIG['MODEL']['model']
    model_name = module_name
    fitfile = model_name+'.fit'
    output_symbolic_func = CONFIG['STAT']['output_symbolic_func']
    
    # Check if fit file exists.
    if not os.path.isfile(fitfile):
        print('ERROR: cannot find the fit file "%s"'%fitfile)
        sys.exit()
        
    # Import model module (in other case, deserialization is not working).
    #model = deserialize_model(model_name+'.model')
    module_path = os.path.join('./',module_name+'.py') # relative path
    module = import_module_by_path(module_name,module_path)
        
    # Get fit object from file.
    f = deserialize_fit(fitfile)
    
    # Print symbolic function.
    if output_symbolic_func:
        print('\n\nSYMBOLIC FUNC (ORIGINAL)=====>')
        print(f.model_final.__symbolic_func__)
        print('\n\nSYMBOLIC FUNC (SUBS)=====>')
        dct = {psym.get_value():f.model_final.__params__[psym.__name__].get_value() \
            for psym in f.model_final.__symbolic_params__} # parameters substitution
        expr_sub = f.model_final.__symbolic_func__.subs(dct) # substitute parameters and inputs 
        print(expr_sub)
            
    # Print fitting statistics.
    print('\nFITTING STAT=====>')
    print('f',f)

    # Print model parameters.
    print('\nMODEL FINAL=====>')
    print(f.__model__)

    # Calculate outlier statistics.
    for grp in f.__fitgroups__:
        if not grp.__active__: continue
        grp.calculate_fit_statistics()

    # Output covariance matrix
    print('\n10-ORDERS OF COVARIANCE ("!" MEANS NEGATIVE VALUE)')
    cc = j.Collection(); 
    for i,row in enumerate(f.__covb__):
        #item = {str(j):e for j,e in enumerate(row)}
        item = {str(j):(str(int(np.log10(abs(e))))+\
            ('' if e>0 else '!')) if j>=i else None for j,e in enumerate(row)} # show only orders of covariance
        item['#'] = i
        cc.update(item)
    order = ['#']+[str(i) for i in range(f.__covb__.shape[0])]
    cc.tabulate(order)

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
            value = np.arange(cmin,cmax,cstep)
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
    resids_weighted = CONFIG['PLOTTING']['resids_weighted']
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

    # X axis for residuals.
    resids_x_axes = CONFIG['PLOTTING']['resids_x_axes']   
    resids_x_axes = resids_x_axes.strip()
    if not resids_x_axes: resids_x_axes=OUTPUT

    # Get model and read parameters from file.
    #model = deserialize_model(modfile)
    model = load_model(CONFIG)
    
    # Get data collection.
    _,col_data = read_fitgroups(CONFIG)
        
    # Do the plotting part.
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
    model_vals = model.calculate(grid)
       
    ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3)
    #plt.title(modfile)
    plt.title('residuals')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False) # set normal format
    leg = []
    plt.plot(resid_axis_data,output_column,'ro'); leg.append('data')
    plt.plot(resid_axis_data,model_vals,'b.'); leg.append('model')

    #if SHOW_N:
    #    for n,v,raxis in zip(nn,output_column,resid_axis_data):
    #        ax1.text(raxis,v,'%d'%n)

    plt.legend(leg)
    plt.grid(True)
    plt.ylabel(OUTPUT)
    
    plt.subplot2grid((4,1),(3,0),rowspan=1,sharex=ax1)
    plt.plot(resid_axis_data,output_column-model_vals,'-')
    plt.xlabel(resids_x_axes)
    plt.legend(['data-model'])
    plt.grid(True)
    
    plt.show()    

COLORS = ['blue','red','green','magenta','black','cyan',
    'yellow','purple','orange','brown','pink','grey']

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
    
    # Get options from the config file.
    #module_name = CONFIG['MODEL']['model']
    #model_name = module_name
    #fitfile = model_name+'.fit'

    # Get argument names.
    argnames = get_arguments(CONFIG)
    n_args = len(argnames)
    
    # Get coordinate bindings.
    bindings = get_argument_bindings(CONFIG)
    bindings_dict = dict(bindings)
    
    # Get gridpes and indexes of unfixed arguments
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

    # Get input and output names
    INPUTS = parse_input_columns(CONFIG)
    OUTPUT = CONFIG['DATA']['output_column']
    
    # Import model module (in other case, deserialization is not working).
    #module_path = os.path.join('./',module_name+'.py') # relative path
    #module = import_module_by_path(module_name,module_path)
    
    # Get fitter object from file.
    #f = deserialize_fit(fitfile)
    
    model = load_model(CONFIG)
    
    # Get fitgroups and calculate statistics.
    fitgroups,_ = read_fitgroups(CONFIG)
    # HOW TO CALCULATE FIT STATISTICS???
    #fitgroups = f.__fitgroups__
    #for grp in fitgroups:
    #    grp.calculate_fit_statistics()
    
    # Get model components list.
    components = get_model_components(CONFIG)
    
    # Do plotting stuff.
    if n_unfixed==2:
        fig = plt.figure()
        #plotter = Axes3D(fig)
        plotter = fig.add_subplot(111, projection='3d')
    elif n_unfixed==1:
        plotter = plt
    else:
        raise NotImplementedError

    print('')

    leg = []

    for fitgroup,grp_color in zip(fitgroups,cycle(COLORS)):
        
        # find index of all data points whose fixed columns correspond to gridspec
        ind = np.full(fitgroup.__inputgrid__.__shape__,True)
        for index_fixed in indexes_fixed:
            argname_fixed,colname_fixed = bindings[index_fixed]
            val_fixed = gridspec_dict[argname_fixed][0]
            ind *= abs(fitgroup[colname_fixed]-val_fixed)<0.001
            
        n_tot = len(ind)
        n_sel = np.count_nonzero(ind)
        print('%s: %d selected out of %d'%(fitgroup,n_sel,n_tot))    
        
        if n_sel==0: continue
        
        leg.append('%s (%d/%d)'%(fitgroup.__name__,n_sel,n_tot))
        
        meshes = fitgroup.__inputgrid__.get_meshes(ind)
        #x_data = meshes[indexes_unfixed[0]]
        #y_data = meshes[indexes_unfixed[1]]
        plot_data = [meshes[i] for i in indexes_unfixed]
        output_data = fitgroup.__output__[ind]
        if plot_outlier_stats:
            stat = fitgroup.stats(outlier_stats_type)[ind]
            stat_title = '(%s)'%outlier_stats_type
            stat_colors = stat
        else:
            stat_title = ''
            #stat_colors = 'black'
            stat_colors = grp_color
        
        if n_unfixed in [1,2]:
            #sc = plotter.scatter(x_data,y_data,output_data,c=stat_colors,
            sc = plotter.scatter(*plot_data,output_data,c=stat_colors,
                        s=marker_size,alpha=scatter_opacity)
        else:
            raise NotImplementedError

    g = Grid(*reduce(lambda x,y:x+y,gridspec))
    meshes = g.get_meshes()
    #x_calc = meshes[indexes_unfixed[0]]
    #y_calc = meshes[indexes_unfixed[1]]
    calc_data = [meshes[i] for i in indexes_unfixed]
    
    if calculate_components:
        results = model.calculate_components(g,compnames=components)
    else:
        results = model.calculate(g); results = [results]
            
    for res,compname,color in zip(results,['MODEL']+components,cycle(COLORS)):
        print('plotting %s (%s)'%(compname,color))
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

    leg += ['calc_model']
    plotter.legend(leg)

    plt.grid(True)
    plt.show()
    
def plot(CONFIG):
    """ interface for the plotting section. Supports a number of standard plots. """
    
    # Get options from the config file.
    plot_mode = CONFIG['PLOTTING']['plot_mode']
    
    # Switch on the plotting options.
    if plot_mode=='residuals':
        plot_residuals(CONFIG)
    elif plot_mode=='sections':
        plot_sections(CONFIG)
    else:
        print('ERROR: unknown plotting mode "%s"'%plot_mode)
        sys.exit()

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
