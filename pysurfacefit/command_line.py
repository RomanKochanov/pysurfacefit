import sys
import argparse

from . import tasks

from .config import print_help, config, template_modules_dict

def set_parameters(config,params):
    params = params.split(';')
    for param in params:
        param = param.strip()
        if not param: continue
        path,val = param.split(':')
        section,par = path.split('.')
        config[section][par] = val

def main():
    """ Main command line driver"""

    parser = argparse.ArgumentParser(description=\
        'Python package for linear and non-linear fitting of multidimensional data.')

    parser.add_argument('--config', type=str, 
        help='Configuration file (mandatory for all steps except --startproject)')

    parser.add_argument('--startproject', type=str, 
        help='Stage 1: start empty project, create dummy config')
        
    parser.add_argument('--template', nargs='+', type=str, 
        help='_________1a: us a template for new project')

    parser.add_argument('--merge', type=str,
        help='_________1b: take defaults from a config file')

    parser.add_argument('--set', type=str,
        help='_________1c: explicitly set project parameters')

    parser.add_argument('--init', dest='init',
        action='store_const', const=True, default=False,
        help='Stage 2: create model package')

    parser.add_argument('--split', dest='split',
        action='store_const', const=True, default=False,
        help='Stage 3: split initial datafile with weighting scheme')

    parser.add_argument('--fit', dest='fit',
        action='store_const', const=True, default=False,
        help='Stage 4: start fitting model to data')

    parser.add_argument('--stat', dest='stat',
        action='store_const', const=True, default=False,
        help='Stage 5: calculate statistics')

    parser.add_argument('--codegen', dest='codegen',
        action='store_const', const=True, default=False,
        help='Stage 6: generate Fortran code for fitted model')

    parser.add_argument('--plot', nargs='*', type=str, 
        help='Stage 7: plot sections of the model(s) and compare to data')

    parser.add_argument('--multicut', nargs='*', type=str, 
        help='Stage 8: plot multiple 1D sections of model(s) and data')

    parser.add_argument('--calc', dest='calc',
        action='store_const', const=True, default=False,
        help='Stage 9: calculate model values on grid')

    args = parser.parse_args() 
    
    CONFIG = config
        
    if args.startproject:
        if not args.template:
            args.template = []
        elif args.template[0]=='help':
            print_help(template_modules_dict)
            sys.exit()       
        for template_pair in args.template:
            module_name,template_name = template_pair.split('.')
            section = getattr(template_modules_dict[module_name],template_name)
            CONFIG.merge_section(section)
        CONFIG['GENERAL']['project'] = args.startproject
        if args.merge:
            CONFIG.load_ini(args.merge,ignore_empty_values=True)
        #if args.set:
        #    set_parameters(CONFIG,args.set)
        # INPUT PARAMETERS
        datafile = input('Enter the name of the data points CSV file: ')
        CONFIG['DATA']['datafile'] = datafile
        input_columns = input('Enter semicolon-separated names for model inputs: ')
        CONFIG['DATA']['input_columns'] = input_columns
        CONFIG['MODEL']['arguments'] = input_columns
        output_column = input('Enter name for model output: ')
        CONFIG['DATA']['output_column'] = output_column
        # Run start project.
        tasks.startproject(CONFIG)
    else:
        if not args.config:
            print('Error: config must be specified (use --config option).')
            sys.exit()
        CONFIG.load_ini(args.config)
        
    if args.set:
        set_parameters(CONFIG,args.set)
    
    if args.init:
        tasks.init(CONFIG)
    elif args.split:
        tasks.split(CONFIG)
    elif args.fit:
        tasks.fit(CONFIG)
    elif args.stat:
        tasks.stat(CONFIG)
    elif args.codegen:
        tasks.codegen(CONFIG)
    elif args.plot is not None:
        if args.plot:
            CONFIG['PLOTTING']['gridspec'] = args.plot[0]
        tasks.plot(CONFIG)
    elif args.multicut is not None:
        if args.multicut:
            CONFIG['MULTICUT']['gridspec'] = args.plot[0]
        tasks.plot_multicut(CONFIG)
    elif args.calc:        
        tasks.calc(CONFIG)
