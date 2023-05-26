import inspect
import numpy as np
import jeanny3 as j
from algopy import UTPM

import copy
from functools import reduce

from ..grids import List

class Model:
    """
    Abstract class for the fit model supporting algorithmic differentiation.
    Description of methods:
        1) self.__func__: the function written in Python, accounts for params and components
        1) self.__calc__: interface between __func__ and the internal machinery
                          __calc__ should be re-defined when making an inherited class (e.g. ModelSympy) 
        3) self.__jac__: usually generates automatically from __func__ or __calc__
    """       
    @property
    def __p__(self): # array of values of active ("unfixed") parameters
        return self.__params__.get_values(active_only=True)
        
    @property
    def __wp__(self): # weights of active ("unfixed") parameters
        return self.__params__.get_weights(active_only=True)
        
    @property
    def __npar__(self):
        return self.__params__.__npar__

    def __call__(self,*inputs):
        """
        Making possible to call the model as a function (magic method).
        Don't confuse with __calc__!!!
        """
        return self.__calc__(self.__params__,*inputs)
        
    def __calc__(self,params,*inputs):
        """
        Interface between the func and Jacobian.
        Don't confuse with __call__!!!
        """
        return self.__func__(params,*inputs)
            
    def calculate(self,grid):
        """ 
        Calculate model on grid defined in pysurface.grids.
        """
        return grid.calculate(lambda *x:self.__calc__(self.__params__,*x))
        
    def calculate_on_grid(self,grid):
        """ 
        Calculate model on grid defined in pysurface.grids.
        """
        return self.calculate(grid)
        
    def calculate_on_collection(self,collection,mapping={}):
        """
        Calculate model on Jeanny3 collection with optional parameter mapping.
        Doesn't support fast calculation for Sympy models.
        Mapping should be list,tuple or dictionary.
        Returns grid and calculated values.
        """
        argspec = inspect.getfullargspec(self.__func__)
        args = argspec.args[2:]
        if type(mapping) in {list,tuple}:
            mapping = {arg:argmap for arg,argmap in zip(args,mapping)}
        map_ = {arg:mapping.get(arg,arg) for arg in args}
        args_ = list(map_.values())
        map_cols = {arg_:c for arg_,c in zip(args_,collection.getcols(args_))}
        grid = List(*reduce(lambda x,y:x+y,[[arg,map_cols[arg]] for arg in map_cols]))
        return grid,self.calculate(grid)

    def assign(self,collection,field,mapping={}):
        """
        Assign values into Jeanny3 collection with optional parameter mapping.
        Supports fast calculation for Sympy models.
        Mapping should be list,tuple or dictionary.
        """
        argspec = inspect.getfullargspec(self.__func__)
        args = argspec.args[2:]
        if type(mapping) in {list,tuple}:
            mapping = {arg:argmap for arg,argmap in zip(args,mapping)}
        map_ = {arg:mapping.get(arg,arg) for arg in args}
        args_ = list(map_.values())        
        cc = collection.getcols(args_+['__ID__'])
        cols = cc[:-1]
        ids = cc[-1]
        map_cols = {arg_:c for arg_,c in zip(args_,cols)}
        grid = List(*reduce(lambda x,y:x+y,[[arg,map_cols[arg]] for arg in map_cols]))
        vals = self.calculate(grid)
        for id_,val in zip(ids,vals):
            collection.getitem(id_)[field] = val

    def calculate_components(self,grid,compnames=[]):
        """
        Calculate model and its components on grid.
        """
        def func_comp(*x):
            res = self.__calc__(self.__params__,*x)
            comps = [self.__components__[name] for name in compnames]
            return [res,*comps]
        hypermesh = grid.calculate(func_comp,dtype=object)
        meshes = [np.empty(hypermesh.shape,dtype=np.float64) for _ in range(len(compnames)+1)]
        for k,_ in enumerate(meshes):
            for i,_ in np.ndenumerate(hypermesh):
                meshes[k][i] = hypermesh[i][k]
        return meshes
        
    def calculate_jac(self,grid):
        """ 
        Calculate model Jacobian (flattened)
        """
        jac_tensor = grid.calculate(lambda *x:self.__jac__(self.__params__,*x),dtype=object,flat=True)
        jac_tensor = np.vstack(jac_tensor)
        return jac_tensor
                   
    def __jac__(self,params,*inputs): 
        # Using algopy direct mode        
        def foo(p):
            params_ = copy.deepcopy(params)
            params_.set_values(p,active_only=True)
            return self.__calc__(params_,*inputs)
        p = UTPM.init_jacobian(list(params.get_values(active_only=True))) 
        foo_val = foo(p)
        algopy_jacobian = UTPM.extract_jacobian(foo_val)
        return algopy_jacobian   
        
    #@abstractmethod
    def __func__(self,params,*inputs):
        raise NotImplementedError
        
    def __units__(self):
        """ 
            Return units of the inputs and output.
            The output must be a nested dict, e.g.:
                {   
                    'inputs': {'r1':'bohr',r2:'bohr',alpha:'deg'},
                    'output': 'cm-1'
                }
        """
        raise NotImplementedError
        
    def save_params(self,filename):
        self.__params__.export_csv(filename)
        
    def load_params(self,filename):
        col = j.import_csv(filename)
        col.index(lambda v:v['name'])
        for ID in col.__dicthash__:
            item = col.__dicthash__[ID]
            if ID in self.__params__.__dicthash__:
                self.__params__.__dicthash__[ID].update(item)
            else:
                self.__params__.__dicthash__[ID] = item
        
    def __repr__(self):
        argspec = inspect.getfullargspec(self.__func__)
        model_args = argspec.args[2:]
        model_name = self.__class__.__name__
        sep = '==================================\n'
        head = 'MODEL CLASS %s: %s\n'%(model_name,', '.join(model_args))
        return sep+head+sep+str(self.__params__)
