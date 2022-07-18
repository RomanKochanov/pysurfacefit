from .base import Model

# ===============================================================================
# ================================ MODEL FORTRAN ================================
# ===============================================================================

# ============
# COMPILERS
# ============

#class Compiler(ABC):
class Compiler:
    """
    Abstract class for compiler.
    """
    def __init__(self,files=[],options=[]):
        self.__files__ = files # wildcards can be used depending on the compiler
        self.__options__ = options # options  to the compiler 
        
    #@abstractmethod
    def run(self):
        raise NotImplementedError
        
    def set_files(self,*files):
        self.__files__ = files
        
    def set_options(self,*options):
        self.__options__ = options
        
class IFort(Compiler):
    """
    Intel Fortran compiler
    """
    def run(self):
        cmd = ' '.join(['ifort',*self.__options__,*self.__files__])
        os.system(cmd) 

# ====================
# SERIALIZABLE DATA        
# ====================
        
def printval(val,type):
    """
    Print the value using the format which 
    minimizes the loss of accuracy.
    """
    if type in {'double precision',
        'float(8)','float(16)','float(32)'}:
        return '%25.15E'%val
    else:
        return str(val)
        
def packvals(array,type,n_per_line=None,by_cols=False,comments=[]):
    """
    Pack values using comma-separated format specifying
    whether to pack by columns (by_cols=True) or by lines (by_cols=False)
    """
    assert len(comments)==0 or len(comments)==len(array),\
        'comments must be either empty or the same length as array'
    if n_per_line is None:
        if type=='integer':
            n_per_line = 20
        elif type in {'double precision','float(8)','float(16)'}:
            n_per_line = 5
        else:
            n_per_line = 3
    if by_cols:
        order = 'F'
    else:
        order = 'C'
    vals = array.flatten(order=order) # flatten using the Fortran order
    buf = ''
    comm = ''
    for i,val in enumerate(vals):
        buf += printval(val,type)
        if i!=len(vals)-1: 
            buf += ', '
        else:
            buf += '  '
        if len(comments)>0: comm += ' '+comments[i]
        #if (i+1)%n_per_line==0 and i!=len(vals)-1:
        if (i+1)%n_per_line==0:
            buf += '  &'
            if len(comments)>0: buf+=' !'+'%3d !'%(i+1)+comm
            buf += '\n'
            comm = ''
    #buf = buf[:-2] # get rid of the last comma and new line
    #buf = buf[:-1] # get rid of the last comma
    return buf
    
class Data:
    """
    Container class for the serializable data items.
    """
    def __init__(self,*items):
        items = list(items)
        self.__items__ = items
        self.__index__ = {item.__name__:i \
            for i,item in enumerate(items)}
        
    def __getitem__(self,name):
        i = self.__index__[item.__name__]
        return self.__items__[i]
        
    def append(self,item):
        if item.__name__ in self.__index__:
            i = self.__index__[item.__name__]
            self.__items__[i] = item
        else:
            self.__index__[item.__name__] = len(self.__items__)
            self.__items__.append(item)
        
    def __get_code__(self):
        code = ''
        for item in self.__items__:
            code += item.__get_code__()
        return code
        
#class DataItem(ABC):
class DataItem:
    """
    Abstract class for the data storage.
    """
    #@abstractmethod
    def __get_code__(self):
        raise NotImplementedError

class FortranValue(DataItem):
    """
    Storage for a scalar Fortran value.
    """
    def __init__(self,name,type,val):
        self.__name__ = name
        self.__fortran_type__ = type
        self.__value__ = val
            
    def __get_code__(self):
        return textwrap.dedent("""
        {type} :: {name} = {value}
        """.format(type=self.__fortran_type__,
                   name=self.__name__,
                   value=self.__value__))
        
class FortranArray(DataItem):
    """
    Storage for a 1D Fortran array
    """
    def __init__(self,name,type,lst,n_per_line=None):
        self.__name__ = name
        self.__fortran_type__ = type
        self.__array__ = np.array(lst)
        self.__n_per_line__ = n_per_line
            
    def __get_code__(self):
        return textwrap.dedent("""
        {type},dimension({n}) :: {name} = &
        """.format(type=self.__fortran_type__,
                   name=self.__name__,
                   n=self.__array__.size)) + \
        '(/&\n' + packvals(self.__array__,self.__fortran_type__,
                        n_per_line=self.__n_per_line__,by_cols=True) +\
        '/)'
        
class FortranMatrix(DataItem):
    """
    Storage for a 2D Fortran array
    """
    def __init__(self,name,type,lst,n_per_line=None):
        self.__name__ = name
        self.__fortran_type__ = type # e.g. "double precision"
        self.__matrix__ = np.array(lst)
        self.__n_per_line__ = n_per_line
    
    def __get_code__(self):
        return textwrap.dedent("""
        {type},dimension({ni}, {nj}) :: {name} = &
        """.format(type=self.__fortran_type__,
                   name=self.__name__,
                   ni=self.__matrix__.shape[0],   # CHECK THE FORTRAN ORDER OF ROWS/COLUMNS!!!
                   nj=self.__matrix__.shape[1])) + \
        '(/&\n' + packvals(self.__matrix__,self.__fortran_type__,
                        n_per_line=self.__n_per_line__,by_cols=True) +\
        '/)'
        
# ===================================================
# TEMPORARY FORTRAN CODE GENERATOR FOR PYTHON MODELS
# ===================================================
# NOTE: __func__ must be manually translated from Python to Fortran!!!
# NOTE2: This class has nothing in common with the rest of the Fortran-oriented 
#        classes defined above!!!

class FortranModule():
    """
    Generate Fortran module from the model object.
    The main function (__func__) is not included since
    it is not trivial to convert Python code to Fortran.
    It attempts to deduce the type and dimensions of the value automatically.
    LIMITATIONS:
        1) Currently, the following data types are supported: 
                int, np.int32, np.int64, 
                float, np.float64, 
                list, np.ndarray
        2) Only one- and two-dimensional arrays are supported.
    """
    def __init__(self,model,item_names):
        # Fitting models
        self.__model__ = model
        # Deduce types and dimensions
        self.__data__ = Data()
        for name in item_names:
            if name[0]=='_': 
                name_ = 'fort'+name # fix leading underscore
                print('WARNING: %s was changed to %s'%(name,name_))
            else:
                name_ = name
            # start deducing...
            val = getattr(model,name)
            if type(val) in {int,np.int32,np.int64}:
                item = FortranValue(name_,'integer',val)
            elif type(val) in {float,np.float64}:
                item = FortranValue(name_,'fouble precision',val)
            elif type(val) in {list,np.ndarray}:
                if type(val)==list: val = np.array(list)
                dtype = val.dtype
                if dtype in {np.dtype('int32'),np.dtype('int64')}: # np.dtype('int32')==np.int32=True, BUT np.dtype('int32') is np.int32=False !!!
                    typestr = 'integer'
                elif dtype in {np.dtype('float32'),np.dtype('float64')}:
                    typestr = 'double precision'
                else:
                    print('dtype>>>',dtype)
                # insert additional item storing the leading dimension of an array/matrix
                item_ = FortranValue('n_'+name_,'integer',val.shape[0])
                self.__data__.append(item_)
                # insert the array/matrix...
                if len(val.shape)==1:
                    item = FortranArray(name_,typestr+'(%d)'%val.shape[0],val)
                elif len(val.shape)==2:
                    item = FortranMatrix(name_,typestr,val.T,n_per_line=val.shape[1])
                else:
                    raise Exception('too many dimensions: ',val.shape)
            self.__data__.append(item)
                
    def generate_code(self):
        # generate code for items
        itemcode = self.__data__.__get_code__()
        
        # generate code for for parameters (__params__)
        params = self.__model__.__params__
        npar = len(params)
        parnames = params.get_names()
        parvalues = params.get_values() 
        parcode = textwrap.dedent("""
        double precision, parameter :: npar = {npar}
        double precision, dimension(npar) :: p = &
        (/&{new_line}{param_values}/)
        """).format(new_line='\n',npar=npar,param_values=\
            packvals(parvalues,type='double precision',n_per_line=1,comments=parnames))
        
        # insert additional useful constants
        constcode = textwrap.dedent("""
        double precision :: PI_IN_RADIANS=3.14159265358979d0
        double precision :: WAVENUMBER_TO_HARTREE=0.4556335d-5
        """)
        
        module_code = textwrap.dedent("""
        module {module_name}
        
        implicit none
        {itemcode}
        {parcode}
        {constcode}        
        end module {module_name}
        """).format(module_name=self.__model__.__class__.__name__,
                itemcode=itemcode,parcode=parcode,constcode=constcode)
        
        return module_code.strip()
        
# =============================
# FORTRAN MODEL IMPLEMENTATION
# =============================
    
class ModelFortran(Model):
    """
    Abstract class for the fit model based on Fortran routines.
    THe model requires Intel Fortran compiler and Adifor (for Jacobian)
    Adding another compilers is planned in the future if this model will be useful.
    The following methods are reimplemented: __func__, __jac__
    The following variables should be set in the __init__ method:
        self.__code__: the function code in a specific format
        self.__data__: supplementary data to be read by the code from the Fortran module
        self.__name__: name, needed to create files on the disc
        self.__params__: python parameter objects to make connection with the Python model
    """
    def __init__(self):
        
        raise NotImplementedError # model fortran is not yet implemeneted
        
        self.__module_f90_name__ = 'module_'+self.__name__
        self.__module_jac_f90_name__ = 'module_'+self.__name__+'__adifor'
        self.__dll_func__ = None # object storing the compiled Fortran module
        self._make_f90_module()
        self._compile_func()
        #self._compile_jac()
    
    def _import(self):
        # import the dll    
        dll_func_name = self.__name__ # in case if BIND(C) option is enabled
        #dll_func_name = 'MODULE_{name}_mp_{name}'.format(name=self.__name__.upper())
        module = ct.CDLL(self.__module_f90_name__+'.dll')
        self.__dll_func__ = getattr(module,dll_func_name)        
        
    def _make_f90_module(self):
        """
        Prepare and save Fortran code for processing by Adifor
        and for importing by the main function.
        """
        # get parameter vector (including fixed ones)
        parvals = self.__params__.get_values()
        
        # split the code into lines 
        func_code = self.__code__
        func_code = textwrap.dedent(func_code)
        func_code = func_code.strip()
        func_code = func_code.split('\n')

        # update the name of the subroutine
        func_code[0] = re.sub('subroutine\s+(.+)\(','subroutine %s('%self.__name__,func_code[0])
        func_code[0] += ' BIND(C, NAME=\'%s\')'%self.__name__  # commenting this will possibly break the code portability

        # insert the export declaration
        func_code.insert(1,'!DEC$ ATTRIBUTES DLLEXPORT :: %s'%self.__name__)

        # get the list of model inputs
        parlist = re.search('subroutine\s+[^\(]+\(([^\(\)]+)\)',func_code[0]).group(1).split(',')
        pars = ','.join(parlist[2:])        
        
        # insert declaration for the interface parameter (assumed-shape arrays give compilation error)
        func_code.insert(2,'double precision, intent(out) :: %s'%parlist[0])
        func_code.insert(2,'double precision, intent(in) :: %s(np__),%s'%(parlist[1],','.join(parlist[2:])))
        
        # get the wrapper code for the function to be used from Fortran
        fortran_wrapper = textwrap.dedent("""
        subroutine wrapper_{name}(res,{pars})
                
        double precision,intent(out) :: res
        double precision,intent(in) :: {pars}
        
        call {name}(res,p,{pars})
        
        end subroutine
        """).format(name=self.__name__,pars=pars)
        
        # join the code lines back
        func_code = '\n'.join(func_code)
        
        # get the complete code of the module
        module_code = textwrap.dedent("""
        module {module_name}
        
        use iso_c_binding
        implicit none
        
        {data}
        {np__}
        {params}
        
        contains
        
        {function_lower_level}
        
        {function_higher_level}
        
        end module {module_name}
        """).format(module_name=self.__module_f90_name__,
            data=self.__data__.__get_code__(),
            np__=FortranValue('np__','integer',len(parvals)).__get_code__(),
            params=FortranArray('p','double precision',
                parvals).__get_code__(),
            function_lower_level=func_code,
            function_higher_level=fortran_wrapper)
        module_code = module_code.strip()
        
        # save the result
        with open(self.__module_f90_name__+'.f90','w') as f:
            f.write(module_code)
        
    def _compile_func(self):
        """
        Compile codes for function and Jacobian, and import the dll modules
        """
        # Compile and import function.
        compiler = IFort(files=[self.__module_f90_name__+'.f90'],options=['/dll','/traceback']) # debug
        #compiler = IFort(files=[self.__module_f90_name__+'.f90'],options=['/dll','/O3']) # release
        compiler.run()
        
    def _compile_jac(self):
        # Compile and import Jacobian.
        compiler = IFort(files=[self.__module_jac_f90_name__+'.f90'],options=['/dll','/O0','/g','/traceback']) # debug
        #compiler = IFort(files=[self.__module_jac_f90_name__+'.f90'],options=['/dll','/O3']) # release
        compiler.run()
          
    #@abstractmethod
    def __func__(self,params,*inputs):
        # make initial preparations
        if self.__dll_func__ is None:
            self._import()
        p = params.get_values()
        # setup the data
        input_pointers = [ct.pointer(ct.c_float(input)) for input in inputs]
        # call the function by passing the ctypes pointer using the numpy function:
        res = ct.c_double(); res_pointer = ct.pointer(res)
        self.__dll_func__(res_pointer,np.ctypeslib.as_ctypes(p),*input_pointers)
        return res.value

    def __jac__(self,params,*inputs):
        raise NotImplementedError()
            
            
# ===============================================================================
# =============================== /MODEL FORTRAN ================================
# ===============================================================================
