import os
import sys
import copy
import subprocess

import numpy as np
import sympy as sy
#import symengine as sy # should be faster than sympy, but not entirely compatible
import numba as nb

from .base import Model

def subprocess_run_crossplatform(command,*args,**kwargs):
    """ Wrapper for subprocess.run to enable basic cross-platform execution.
    ==============
    sys.platform:
    ==============
    Linux 	'linux'
    Windows 	'win32'
    Windows/Cygwin 	'cygwin'
    Mac OS X 	'darwin'
    """
    if type(command) not in (list,tuple):
        command = [command]
    exefile = command[0]
    if sys.platform in ['linux','darwin','cygwin']:
        exefile = './'+exefile
    subprocess.run([exefile]+command[1:],*args,**kwargs)
        
# ===============================================================================
# ================================ MODEL SYMPY ==================================
# ===============================================================================

# HOW TO GET THE VAR NAMES OF THE FUNCTION:
# https://stackoverflow.com/questions/582056/getting-list-of-parameter-names-inside-python-function

def get_argnames(foo):
    argcount = foo.__code__.co_argcount
    argnames = foo.__code__.co_varnames[:argcount]
    return argnames

def gradient(scalar_function,variables):
    """https://stackoverflow.com/questions/21166771/is-there-a-vectorized-way-to-calculate-the-gradient-in-sympy?rq=1"""
    matrix_scalar_function = sy.Matrix([scalar_function])
    return matrix_scalar_function.jacobian(sy.Matrix([variables]))    

def create_pots_interface():
    """ 
        Create POTS interface for the DVR program only.
        POTS is calling PES. Generally, this should not depend on the model
        and custom parameter usage, such as nuclear masses.
    """
    
    pots_source = \
    'subroutine pots(vout,r1,r2,alpha)' + '\n' + \
    '' + '\n' + \
    'implicit none' + '\n' + \
    '' + '\n' + \
    '! ====================================================================================' + '\n' + \
    '! POTS is the interface between the DVR3D code and the user defined PES routine PES' + '\n' + \
    '!' + '\n' + \
    '! All data are transferred into POTS either through the argument list (r1,r2,alpha) or' + '\n' + \
    '! common blocks' + '\n' + \
    '!' + '\n' + \
    '! ====================================================================================' + '\n' + \
    '' + '\n' + \
    '! This particular case is the interface with the PS mass dependent water PES ' + '\n' + \
    '' + '\n' + \
    'double precision, intent(in) :: r1    ! first bondlength [BOHR]' + '\n' + \
    'double precision, intent(in) :: r2    ! second bondlength [BOHR]' + '\n' + \
    'double precision, intent(in) :: alpha ! bondangle [RADIAN]' + '\n' + \
    'double precision, intent(out) :: vout    ! returned PES value [AU]' + '\n' + \
    '' + '\n' + \
    '! ==================================================================================' + '\n' + \
    '! ==================================================================================' + '\n' + \
    '' + '\n' + \
    'common /mass/ xmass(3) ! vib. masses (must be in in a.u.)' + '\n' + \
    'real*8 xmass' + '\n' + \
    '' + '\n' + \
    'double precision pes' + '\n' + \
    '' + '\n' + \
    '! call PES function to calculate value of PES in' + '\n' + \
    '! the geometry point {r1,r2,alpha} ' + '\n' + \
    '' + '\n' + \
    'vout = pes(.TRUE.,0,r1,r2,alpha,xmass)' + '\n' + \
    '' + '\n' + \
    'end'
    
    with open('pots.f90','w') as fil:
        fil.write(pots_source)

def create_pes_interface(model):
    """ 
        Create an interface for the VTET/DVR programs     
        Source file must correspond to the models units.
        
        Note that, depending on the "au" parameters,
        the units of the input arguments (and output value),
        supplied in the "pes" interface, are changing:
            if au=true: r supplied in bohr, pes expected in hatree
            if au=false: r supplied in angstroem, pes expected in wavenumbers
        
        In both cases alpha supplied in radians.
    """
    units = model.__units__()
    input_units = units['input']
    output_units = units['output']
    input_names = model.__input_names__
    
    def convert_units_pes_to_model(quantity,model_unit):
        """
        All parameters are strings.
        quantity - paramter to convert units for
        model_units - units of the input in Python model __func__
        """
        if model_unit.lower() in ['bohr','au','a.u.']:
            return quantity
        elif model_unit.lower() in ['ang','ang.','a','angstroem','angstrem']:
            return quantity+'*BOHR_TO_ANGSTROEM'
        elif model_unit.lower() in ['rad','rad.','radian','radians']:
            return quantity
        elif model_unit.lower() in ['deg','deg.','degree','degrees']:
            return quantity+'/DEGREE_TO_RADIAN'
        elif model_unit.lower() in ['hartree']:
            return quantity+'/WAVENUMBER_TO_HARTREE'
        elif model_unit.lower() in ['cm-1','wavenumber','wavenumbers']:
            return quantity
        else:
            raise Exception('unknown units for quantity "%s": %s'%(quantity,model_unit))
    
    pes_source = \
    'function pes(au,ido,r1,r2,r3,mass)' + '\n\n' + \
    '' + \
    'implicit none' + '\n\n' + \
    '' + \
    'logical au' + '\n' + \
    'integer ido' + '\n' + \
    'real(8) pes,r1,r2,r3,mass(*)' + '\n' + \
    'real(8) r1_,r2_,r3_' + '\n\n' + \
    '' + \
    'real(8) mass1,mass2' + '\n' + \
    'real(8) BOHR_TO_ANGSTROEM,PI_IN_RADIANS,DEGREE_TO_RADIAN,WAVENUMBER_TO_HARTREE' + '\n\n' + \
    '' + \
    'intent(in) :: au,ido,r1,r2,r3' + '\n\n' + \
    '' + \
    'parameter(BOHR_TO_ANGSTROEM=0.529177249d0)' + '\n' + \
    'parameter(PI_IN_RADIANS=3.14159265358979d0)' + '\n' + \
    'parameter(DEGREE_TO_RADIAN=1.74532925199433d-2)' + '\n' + \
    'parameter(WAVENUMBER_TO_HARTREE=0.4556335d-5)' + '\n\n' + \
    '' + \
    '!===========================================================' + '\n' + \
    '! Interface routine between VTET and PES codes' + '\n' + \
    '!===========================================================' + '\n' + \
    '' + \
    'real(8) pes1' + '\n\n' + \
    '' + \
    '! UNITS1: Convert input units according to the old pes_ expectations.' + '\n' + \
    '! au==TRUE => r in bohr, theta in rad' + '\n' + \
    '! au==FALSE => r in angstroem, rho in rad' + '\n' + \
    'if (au) then' + '\n' + \
    '    r1_ = r1' + '\n' + \
    '    r2_ = r1' + '\n' + \
    '    r3_ = r3' + '\n' + \
    'else' + '\n' + \
    '    r1_ = r1/BOHR_TO_ANGSTROEM' + '\n' + \
    '    r2_ = r1/BOHR_TO_ANGSTROEM' + '\n' + \
    '    r3_ = PI_IN_RADIANS-r3' + '\n' + \
    'endif' + '\n\n' + \
    '' + \
    '! UNITS2: Convert input units according to the generated Python model.' + '\n' + \
    'r1_ = %s'%convert_units_pes_to_model('r1_',input_units[input_names[0]]) + '\n' + \
    'r2_ = %s'%convert_units_pes_to_model('r2_',input_units[input_names[1]]) + '\n' + \
    'r3_ = %s'%convert_units_pes_to_model('r3_',input_units[input_names[2]]) + '\n\n' + \
    '' + \
    'mass1 = mass(1)' + '\n' + \
    'mass2 = mass(2)' + '\n\n' + \
    '' + \
    'call sympy_generated(pes1,r1_,r2_,r3_%s)'%(',mass1,mass2' if len(input_names)==5 else '') + '\n\n' + \
    '' + \
    '! UNITS3: Convert output units according to the old pes_ expectations.' + '\n' + \
    'pes1 = %s'%convert_units_pes_to_model('pes1',output_units) + '\n\n' + \
    '' + \
    '! UNITS4: Convert input units according to the generated Python model.' + '\n' + \
    'if (au) then' + '\n' + \
    '    pes1 = pes1*WAVENUMBER_TO_HARTREE' + '\n' + \
    'endif' + '\n\n' + \
    '' + \
    'if(ido==0) then' + '\n' + \
    '    pes=pes1' + '\n' + \
    'else' + '\n' + \
    '    stop "ERROR: noadifor version does not imply PES differentiation"' + '\n' + \
    'endif' + '\n\n' + \
    '' + \
    'end'
    
    with open('pes.f90','w') as fil:
        fil.write(pes_source)

def compile_gridcalc_fortran(model,compiler='gfortran'):
    """
    Compile the gridcalc with a model given in filename.
    The compiler executable must be accessible!
    """

    filename = model.__f90_filename__
    params_list = ','.join(model.__input_names__)
    exefile = model.__class__.__name__

    # Define the gridcalc function to return.
    def gridcalc(exefile,grid):
        # create grid file
        meshes = grid.get_meshes(flat=True)
        #meshes = [mesh.flatten() for mesh in meshes]
        gridfile = '~.pts'
        #np.savetxt(gridfile,list(zip(r1_,r2_,alpha_)))
        np.savetxt(gridfile,list(zip(*meshes)))
        # call to gridcalc program
        #subprocess.run(['./'+exefile,gridfile])
        subprocess_run_crossplatform([exefile,gridfile])
        outputfile = gridfile+'.out'
        data = np.loadtxt(outputfile)
        data = list(zip(*data))
        res = np.reshape(data[-1],grid.__shape__)
        return res
        
    # Create the gridcalc main source file.
    gridcalc_main_source = \
        'program gridcalc' + '\n\n' + \
        '!use ifport' + '\n\n' + \
        'implicit none' + '\n\n' + \
        'character gridfilename*100' + '\n' + \
        'integer ios' + '\n' + \
        'integer ch_in,ch_out' + '\n' + \
        'double precision :: %s,res'%params_list + '\n\n' + \
        'call getarg(1,gridfilename)' + '\n\n' + \
        'print *,\'gridfilename \',gridfilename' + '\n\n' + \
        'open(unit=1,file=gridfilename,action=\'READ\',status=\'OLD\')' + '\n' + \
        'open(unit=2,file=trim(gridfilename)//\'.out\',action=\'WRITE\',status=\'REPLACE\')' + '\n\n' + \
        'print *,\'start loop...\'' + '\n\n' + \
        'do while (.true.)' + '\n' + \
        '  read(1,*,iostat=ios) %s'%params_list + '\n' + \
        '  if( ios > 0 ) then' + '\n' + \
        '    stop "problem somewhere"' + '\n' + \
        '  else if( ios < 0 ) then' + '\n' + \
        '    exit' + '\n' + \
        '  endif' + '\n' + \
        '  call sympy_generated(res,%s)'%params_list + '\n' + \
        '  write(2,\'(4E25.15)\') %s,res'%params_list + '\n' + \
        'end do' + '\n' + \
        '200 continue' + '\n' + \
        'print *,\'end loop\'' + '\n\n' + \
        'close(1)' + '\n' + \
        'close(2)' + '\n\n' + \
        'end program gridcalc'
        
    with open('gridcalc.f90','w') as fil:
        fil.write(gridcalc_main_source)  
        
    # Compile the existing Fortran source.
    subprocess.run([compiler,'-o',exefile,filename,'gridcalc.f90'])
    #subprocess_run_crossplatform([compiler,'-o',exefile,filename,'gridcalc.f90'])
    
    # Get exefile and create the resulting function.
    filestem,_ = os.path.splitext(filename)
    foo = lambda grid: gridcalc(exefile,grid)
    
    return foo
    
        
# ATTENTION =================================
# TOO MANY PARAMETERS WILL SPAWN ERROR IN LAMBDIFY:
# SyntaxError: more than 255 arguments
# POSSIBLE SOLUTION IS HERE:
# https://stackoverflow.com/questions/45220567/sympy-lambdify-function-with-large-array-input    
    
class ModelSympy(Model):
    """
    Abstract class for the fit model based on Sympy symbolic library.
    Made for two major reasons, mostly for being able to:
        1) ...compute Jacobian efficiently
        2) ...export the model and its parameters to the Fortran module
    The following variables should be set in the __init__ method:
        self.__params__: python parameter objects to make connection with the Python model
    """
    ##########################################################
    ################# CALCULATE FUNCTION #####################
    ##########################################################
    
    def __sympy_initialize_func__(self):
        print('')
        print('==========================================')
        print('<<<< calling __sympy_initialize_func__ >>>>')
        print('==========================================')
        print('Progress:')
        
        # Create Sympy objects for inputs.
        print('     - creating sympy objects for inputs')
        argnames = get_argnames(self.__func__); argnames = list(argnames)
        argnames = argnames[2:] # omit "self" and "params" arguments
        inputs_ = [sy.Symbol(argname,real=True) for argname in argnames]
        self.__symbolic_inputs__ = inputs_
        self.__input_names__ = argnames
        
        # Create Sympy objects for parameters.
        print('     - creating sympy objects for parameters')
        #sympy_pars = [sy.Symbol(p.__name__,real=True) for p in self.__params__ if p.__flag__] # ONLY ACTIVE PARAMETERS
        sympy_pars = [sy.Symbol(p.__name__,real=True) for p in self.__params__]
        params_ = copy.deepcopy(self.__params__)
        #params_.set_values(sympy_pars,active_only=True) # ONLY ACTIVE PARAMETERS
        params_.set_values(sympy_pars,active_only=False)
        self.__symbolic_params__ = params_
        
        # Get the Sympy expression by calling function with the Sympy objects.
        print('     - get the Sympy expression by calling function with the Sympy objects')
        self.__symbolic_func__ = self.__func__(params_,*inputs_)
        
        # Create lambdified Python function from sympy expression.
        print('     - create lambdified Python function from sympy expression')
        args_lambdify = [psym.get_value() for psym in self.__symbolic_params__]+self.__symbolic_inputs__
        self.__lambdified_func__ = sy.lambdify(args_lambdify,self.__symbolic_func__)
        
        # Create compiled (numbified) code from the lambdified function.
        print('     - create compiled (numbified) code from the lambdified function')
        self.__numbified_func__ = nb.njit(self.__lambdified_func__)
        
        print('==========================================\n')
                
    def __calc_symbolic__(self,params,*inputs): # symbolic expression->numeric result (through substitution) ONLY FOR DEBUGING PURPOSES, VERY SLOW!!!
        dct = {psym.get_value():params[psym.__name__].get_value() for psym in self.__symbolic_params__} # parameters substitution
        dct.update({v:inp for v,inp in zip(self.__symbolic_inputs__,inputs)})
        expr_sub = self.__symbolic_func__.subs(dct) # substitute parameters and inputs        
        return expr_sub.evalf()

    def __calc_lambdified__(self,params,*inputs): # lambdified expression->numeric result
        args = [p.get_value() for p in params]+list(inputs)
        return self.__lambdified_func__(*args)

    def __calc_numbified__(self,params,*inputs): # compiled lambdified expression->numeric result
        args = [p.get_value() for p in params]+list(inputs)
        return self.__numbified_func__(*args)
        
    def calculate_components(self,grid,compnames=[]): # calculate model and it's components on grid
        def func_comp(*x):
            # get the symbolic components
            res = self.__func__(self.__params__,*x)
            comps = [self.__components__[name] for name in compnames]
            # evaluate all components
            #res = res.evalf()   # WORKS FINE WITHOUT THIS??????
            #comps = [comp.evalf() for comp in comps]   # WORKS FINE WITHOUT THIS??????
            # return evaluated components
            return [res,*comps]
        hypermesh = grid.calculate(func_comp,dtype=object)
        meshes = [np.empty(hypermesh.shape,dtype=np.float64) for _ in range(len(compnames)+1)]
        for k,_ in enumerate(meshes):
            for i,_ in np.ndenumerate(hypermesh):
                meshes[k][i] = hypermesh[i][k]
        return meshes
        
    def __calc__(self,params,*inputs):
        # check for the model initialization
        if '__initialized_func__' not in self.__dict__ or self.__initialized_func__ is False:
            self.__sympy_initialize_func__()
            self.__initialized_func__ = True
            
        # check the validity of algorithm on each N_th calculation
        if '__calc_counter__' not in self.__dict__:
            self.__calc_counter__ = 0
        else:
            self.__calc_counter__ += 1
        N_CHECK = 100000 # have it sufficiently large to avoid littering the output
        #DO_CHECK = True
        if '__check_symbolic__' not in self.__dict__:
            DO_CHECK = True
        else:
            DO_CHECK = self.__check_symbolic__
        if DO_CHECK and self.__calc_counter__%N_CHECK==0:
            print('CALL_COUNTER: %d'%self.__calc_counter__)
            self.__check__(params,*inputs)
        
        # call on of the calc functions depending on an internal switch
        if self.__calc_switch__ == 'symbolic':
            return self.__calc_symbolic__(params,*inputs)
        elif self.__calc_switch__ == 'lambdified':
            return self.__calc_lambdified__(params,*inputs)
        elif self.__calc_switch__ == 'numbified':
            return self.__calc_numbified__(params,*inputs)
        else:
            raise Exception('unknown calc option "%s"'%self.__calc_switch__)

    ##########################################################
    ################# CALCULATE JACOBIAN #####################
    ##########################################################
        
    def __sympy_initialize_jac__(self):
        print('')
        print('==========================================')
        print('<<<< calling __sympy_initialize_jac__ >>>>')
        print('==========================================')
        print('Progress:')
                
        # Get the Sympy expression for jacobian from __symbolic_func__.
        print('     - get the Sympy expression by calling function with the Sympy objects')
        #sympy_pars = [spar.get_value() for spar in self.__symbolic_params__] # ALL PARAMETERS (CAN BE INEFFICIENT IF THE ACTIVE SUBSET IS MUCH SMALLER THAT THE TOTAL PARAMETER SET)
        sympy_pars = [spar.get_value() for spar in self.__symbolic_params__ if spar.__flag__] # !!! ONLY ACTIVE PARAMETERS (CAN CAUSE ERRORS IF THE ACTIVE SUBSET CHANGES WITHIN ONE SESSION)
        self.__symbolic_jac__ = gradient(self.__symbolic_func__,sympy_pars) # symbolic jacobian with respect to all parameters
        
        # Create lambdified Python function from sympy expression.
        print('     - create lambdified Python function from sympy expression')
        args_lambdify = [psym.get_value() for psym in self.__symbolic_params__]+self.__symbolic_inputs__
        self.__lambdified_jac__ = sy.lambdify(args_lambdify,self.__symbolic_jac__)
        
        
        # Create compiled (numbified) code from the lambdified function.
        print('     - create compiled (numbified) code from the lambdified function')
        self.__numbified_jac__ = nb.njit(self.__lambdified_jac__)
        
        print('==========================================\n')
        
    def __jac_symbolic__(self,params,*inputs): # symbolic expression->numeric result (through substitution) ONLY FOR DEBUGING PURPOSES, VERY SLOW!!!
        dct = {psym.get_value():params[psym.__name__].get_value() for psym in self.__symbolic_params__} # parameters substitution
        dct.update({v:inp for v,inp in zip(self.__symbolic_inputs__,inputs)})
        expr_sub = self.__symbolic_jac__.subs(dct) # substitute parameters and inputs  
        return np.array(expr_sub.tolist()[0],dtype=np.float64) # using tolist() method; is there a better way to do this?

    def __jac_lambdified__(self,params,*inputs): # lambdified expression->numeric result
        args = [p.get_value() for p in params]+list(inputs)
        return self.__lambdified_jac__(*args)[0]

    def __jac_numbified__(self,params,*inputs): # compiled lambdified expression->numeric result
        args = [p.get_value() for p in params]+list(inputs)
        return self.__numbified_jac__(*args)[0]
        
    def __jac__(self,params,*inputs): 
        # check for the jacobian initialization
        if '__initialized_jac__' not in self.__dict__ or self.__initialized_jac__ is False:
            self.__sympy_initialize_jac__()
            self.__initialized_jac__ = True
            
        # call on of the calc functions depending on an internal switch
        if self.__calc_switch__ == 'symbolic':
            return self.__jac_symbolic__(params,*inputs)
        elif self.__calc_switch__ == 'lambdified':
            return self.__jac_lambdified__(params,*inputs)
        elif self.__calc_switch__ == 'numbified':
            return self.__jac_numbified__(params,*inputs)
        else:
            raise Exception('unknown calc option "%s"'%self.__calc_switch__)            
            
    ##########################################################
    ############# CHECK FUNCTION CALCULATION  ################
    ##########################################################
        
    def __check__(self,params,*inputs):
        """
        Check the approximate equality of results for four methods of
        the function computations:
        1) Via __symbolic_func__: pure symbolic calculation though substitution
        2) Via __lambdified_func__: pure python function derived from __symbolic_func__
        3) Via __numbified_func__: __symbolic_func__ compiled with Numba
        4) Via __func__: initial function (using numerical parameters and inputs and symbolic functions)
        """
        print('')
        print('====================================================')
        print('<<<< checking the validity of symbolic approach >>>>')
        print('====================================================')
        print(' 1) Via __symbolic_func__: pure symbolic calculation though substitution')
        print(' 2) Via __lambdified_func__: pure python function derived from __symbolic_func__')
        print(' 3) Via __numbified_func__: __lambdified_func__ compiled with Numba')
        print(' 4) Via __func__: initial function (using numerical parameters and inputs and symbolic functions)')
        print('')
        res__symbolic_func__ = self.__calc_symbolic__(params,*inputs)
        res__lambdified_func__ = self.__calc_lambdified__(params,*inputs)
        res__numbified_func__ = self.__calc_numbified__(params,*inputs)
        res__func__ = self.__func__(params,*inputs)
        vals = [res__symbolic_func__, res__lambdified_func__, res__numbified_func__, res__func__]
        n = len(vals)
        diffs = np.zeros([n,n])
        abs_errors = []
        for i in range(n):
            val_i = vals[i]
            for j in range(i,n):
                val_j = vals[j]
                diff = val_i - val_j
                diffs[i,j] = diff
                abs_errors.append(np.abs(diff))
        min_ae = min(abs_errors); max_ae = max(abs_errors)
        print('MAX_ABS_ERROR: %f; MIN_ABS_ERROR: %f\n'%(min_ae,max_ae))
        print('RESIDUALS (i-j):')
        print(tab([[i]+list(line) for i,line in enumerate(diffs)],headers=range(n+1)))
        #THRESH = 1.0E-10
        THRESH = 1.0E-10
        QUIT_ON_ERROR = False
        if max_ae>THRESH:
            msg = 'failed with absolute error maximum of %e (greater than %e)'%(max_ae,THRESH)
            if QUIT_ON_ERROR:                
                raise Exception(msg)
            else:
                warn(msg)
        else:
            print('...PASSED VALIDITY TEST')
        print('====================================================\n')

    ##########################################################
    ############### FORTRAN CODE GENERATION  #################
    ##########################################################
            
    def generate_fortran_90(self):    # BACKUP
        """Generate Fortran 90 code from the model with current parameters"""
        subroutine_name = 'sympy_generated'
        result_var = 'calc__result__'
        
        # Create the Sympy expression
        dct = {psym.get_value():self.__params__[psym.__name__].get_value() for psym in self.__symbolic_params__} # parameters substitution
        expr_sub = self.__symbolic_func__.subs(dct) # substitute parameters and inputs        
        expr_sub_f90 = sy.fcode(expr_sub.evalf(),assign_to=result_var,standard=90,source_format='free')  # get rid of functions of constants;  

        # Create the raw Fortran code
        code = \
        'subroutine %s(%s)\n\n'%(subroutine_name,','.join([result_var]+self.__input_names__))+\
        ''.join(['double precision :: %s\n'%inpname for inpname in self.__input_names__])+\
        'double precision :: %s\n\n'%result_var+\
        expr_sub_f90+'\n\n'+\
        'end subroutine %s'%subroutine_name
        
        # Save the code to file.
        f90_filename = self.__class__.__name__+'.f90'
        print('\n==========================')
        print('CREATING THE FORTRAN CODE')
        print('==========================')
        print('   saving to file "%s"'%f90_filename)
        with open(f90_filename,'w') as fout:
            fout.write(code)
        
        # Save filename
        self.__f90_filename__ = f90_filename
        
        # Compile the grid function
        #self.__function_fortran__ = compile_gridcalc_fortran(self,f90_filename)
        
    def compile_fortran_90(self,compiler='gfortran'):
        self.__function_fortran__ = compile_gridcalc_fortran(self,compiler)
        
        #return code

#    def generate_fortran_90_with_vars(self): # CAN BE DELETED?
#        """Generate Fortran 90 code from the model with current parameters WITH ADDITIONAL VARIABLE DECLARATIONS"""
#        #subroutine_name = self.__class__.__name__
#        subroutine_name = 'sympy_generated'
#        #result_var = subroutine_name+'__result__'
#        result_var = 'calc__result__'
#        expr = self.__symbolic_func__ # original expression        
#        #expr_f90 = sy.fcode(expr,assign_to=result_var,standard=90,source_format='free')        
#        expr_f90 = sy.fcode(expr.evalf(),assign_to=result_var,standard=90,source_format='free')   # get rid of functions of constants; 
#        code = \
#        'subroutine %s(%s)\n\n'%(subroutine_name,','.join([result_var]+self.__input_names__))+\
#        ''.join(['double precision :: %s\n'%inpname for inpname in self.__input_names__])+\
#        'double precision :: %s\n\n'%result_var+\
#        ''.join(['double precision :: %s=%s\n'%\
#            (par.__name__,('%25.17E'%par.__value__).replace('E','D')) for par in self.__params__])+\
#        '\n'+expr_f90+'\n\n'+\
#        'end subroutine %s'%subroutine_name
#        self.__function_fortran__ = compile_gridcalc_fortran(self,expr_sub_f90)
#        raise NotImplementedError # create a stub for this method
#        #return code
        
# ===============================================================================
# ================================ MODEL SYMPY ==================================
# ===============================================================================
