import numpy as np
from tabulate import tabulate as tab

class Parameters:
    """
    Aggregator for fitting parameters.
    This is not necessary to use in models: one can just use the arrays for parameter.
    """
    def __init__(self,N=0,prefix=None,names=None,
                 values=None,flags=True,mins=None,maxs=None,
                 weights=None,group='default'):
        self.__objects__ = []
        self.__index__ = {}
        self.generate(N=N,prefix=prefix,group=group,
                      names=names,values=values,flags=flags,
                      mins=mins,maxs=maxs,weights=weights)
                      
    def dump(self): # (DE)SERIALIZATION FOR PAPER
        parlist = []
        for par in self.__objects__:
            parlist.append(par.dump())
        return parlist
        
    def load(self,parlist): # (DE)SERIALIZATION FOR PAPER
        self.__objects__ = []
        self.__index__ = {}
        for dct in parlist:
            par = Par(name=dct['name']); par.load(dct)
            self.append(par)        
        
    def getcode(self,indent=4,attrs=None): # (DE)SERIALIZATION FOR PAPER
        buf = ''
        for par in self.__objects__:
            buf += (' '*indent+'%s,\n')%par.getcode(attrs=attrs)
        return '[\n%s]'%buf
        
    def __iter__(self):
        return iter(self.__objects__) # strange peculiarity of Python 3
        
    def __len__(self):
        return len(self.__objects__)
                      
    @property
    def __npar__(self):
        return len(self.__objects__)
        
    def __getitem__(self,name):
        """
        Get parameter object by its name.
        """
        return self.__objects__[self.__index__[name]]
        
    def append(self,pars,group='default'):
        """
        Append a parameter to the set. Pars can be either a parameter,
        or a list of parameters.
        """
        if type(pars) not in {list,tuple,Parameters}:
            pars = [pars]
        for par in pars:
            par.__group__ = group
            if par.__name__ not in self.__index__:
                self.__index__[par.__name__] = len(self.__objects__)
                self.__objects__.append(par)
            else:
                self.__objects__[self.__index__[par.__name__]] = par
                        
    def generate(self,N=None,prefix=None,names=None,group='default',
                 values=None,flags=True,mins=None,maxs=None,weights=None):
        """
        Generate n parameters with predefined values.
        """
        if N is None and names is None:
            raise Exception('both N and names are None')
        if prefix is None: prefix = ''
        if names is None:
            # generate names:
            names = [prefix+'%d'%i for i in range(N)]
        else:
            if len(names)!=N:
                raise Exception('len(names)!=N')
        N = len(names)
        if type(weights) not in {list,tuple}:
            weights = [weights for i in range(N)]
        if type(values) not in {list,tuple}:
            values = [values for i in range(N)]
        if type(flags) not in {list,tuple}:
            flags = [flags for i in range(N)]
        if type(mins) not in {list,tuple}:
            mins = [mins for i in range(N)]
        if type(maxs) not in {list,tuple}:
            maxs = [maxs for i in range(N)]
        for name,weight,value,flag,min,max in \
                zip(names,weights,values,flags,mins,maxs):
            self.append(Par(name=name,value=value,flag=flag,
                    min=min,max=max,weight=weight),
                    group=group)

    def get_groups(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array([self.__objects__[i].__group__ for i in ii])
                    
    def get_names(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array([self.__objects__[i].__name__ for i in ii])
        
    def get_values(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array([self.__objects__[i].__value__ for i in ii])

    def set_values(self,values,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)        
        for i,val in zip(ii,values):
            self.__objects__[i].__value__ = val
        
    def get_bounds(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return [self.__objects__[i].get_bounds() for i in ii]
        
    def get_weights(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array([self.__objects__[i].__weight__ for i in ii])
    
    def get_flags(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array([self.__objects__[i].__flag__ for i in ii])
        
    def get_nonzero_weights_and_values(self):
        # THIS CAN HELP FITTING BY REMOVING THE ZERO-WEIGHTED COMPONENTS IN THE RUBBER
        # TO BE IMPLEMENTED
        raise NotImplementedError
        
    def get_index(self,active_only=False,group=None):
        # REDO!!!
        if active_only and group:
            ii = [i for i,p in enumerate(self.__objects__) if p.__group__==group and p.__flag__]
        elif active_only and not group: 
            ii = [i for i,p in enumerate(self.__objects__) if p.__flag__]  
        elif not active_only and group:
            ii = [i for i,p in enumerate(self.__objects__) if p.__group__==group]            
        elif not active_only and not group:
            ii = [i for i,p in enumerate(self.__objects__)]
        return ii
        
    def __repr__(self):
        groups = self.get_groups()
        names = self.get_names()
        values = self.get_values()
        bounds = self.get_bounds()
        weights = self.get_weights()
        flags = self.get_flags()
        indices = []; i = 0
        for flag in flags:
            if flag:
                indices.append(i); i += 1
            else:
                indices.append('')
        return tab(zip(indices,groups,names,values,bounds,weights,flags),
                headers=('#','group','names','values','bounds','weights','flags'))

def dumpval(val):
    """ Simple alternative to json.dumps for VERY limited set of value types:
        1) only scalar types are supported: int, float, str, bool
    """
    if type(val) is str:
        return '"%s"'%val
    else:
        return str(val)
                
class Par:
    """
    Fitting parameter storing information about initial values, bounds, rubber weight, and active flag.
    Each parameter has its unique name.
    """
    MAP = { # (DE)SERIALIZATION FOR PAPER
        'name': '__name__',
        'group': '__group__',
        'value': '__value__',
        'flag': '__flag__',
        'min': '__min__',
        'max': '__max__',
        'weight': '__weight__',
    }
    
    def __init__(self,name,value=None,flag=None,min=None,max=None,weight=None,group=None):
        self.__name__ = name
        self.__group__ = 'default' if group is None else group
        self.__value__ = np.random.random() if value is None else value
        self.__flag__ = False if flag is None else flag
        self.__min__ = -np.inf if min is None else min
        self.__max__ = np.inf if max is None else max
        self.__weight__ = 0.0 if weight is None else weight
        
    def dump(self): # (DE)SERIALIZATION FOR PAPER
        keys = set(self.MAP)
        return {name:getattr(self,self.MAP[name]) for name in keys}
        
    def load(self,dct): # (DE)SERIALIZATION FOR PAPER
        for name in dct:
            if name in self.MAP: setattr(self,self.MAP[name],dct[name])
                
    def getcode(self,attrs=None): # LIMITED AND DISPOSABLE ALTERNATIVE OF THE PROPER DUMP FUNCTIONAL; SHOULD BE USED FOR THE OZONE PAPER ONLY!!!
        if attrs is None: attrs = list(self.MAP.keys())
        dct = {name:getattr(self,self.MAP[name]) for name in attrs}
        return '{%s}'%(', '.join(['"%s":%s'%(name,dumpval(dct[name])) for name in dct]))
                              
    def get_value(self):
        return self.__value__
        
    def set_value(self,value):
        self.__value__ = value
        
    def get_flag(self):
        return self.__flag__
        
    def set_flag(self,flag):
        self.__flag__ = flag        
        
    def get_bounds(self):
        return [self.__min__,self.__max__]
        
    def set_bounds(self,bnds):
        self.__min__ = bnds[0]
        self.__max__ = bnds[1]
        
    def get_min(self):
        return self.__min__
        
    def set_min(self,min):
        self.__min__ = min

    def get_max(self):
        return self.__max__
        
    def set_max(self,max):
        self.__max__ = max
        
    def get_weight(self):
        return self.__weight__
        
    def set_weight(self,wht):
        self.__weight__ = wht
        
    def get_name(self):
        return self.__name__    
        
    def __repr__(self):
        return '{}:{}:{}({}/[{}:{}])*{}'.format(self.__group__,
        self.__name__,self.__flag__,self.__value__,self.__min__,
        self.__max__,self.__weight__)
