import numpy as np

from algopy import UTPM

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
            
    def calculate(self,grid): # calculate model on grid defined in grids.py
        return grid.calculate(lambda *x:self.__calc__(self.__params__,*x))

    def calculate_components(self,grid,compnames=[]): # calculate model and it's components on grid
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
        
    def calculate_jac(self,grid):# returns flattened array
        jac_tensor = grid.calculate(lambda *x:self.__jac__(self.__params__,*x),dtype=object,flat=True)
        jac_tensor = np.vstack(jac_tensor)
        return jac_tensor
                   
    def __jac__(self,params,*inputs): 
        # Using algopy direct mode        
        def foo(p):
            params_ = copy.deepcopy(params)
            params_.set_values(p,active_only=True)
            return self.__calc__(params_,*inputs)
        t = time()
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
        
    def __repr__(self):
        return str(self.__params__)
