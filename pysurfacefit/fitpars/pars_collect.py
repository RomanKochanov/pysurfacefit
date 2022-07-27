""" New version of the fit Parameters based on Jeanny collections."""

import numpy as np
from tabulate import tabulate as tab

from jeanny3 import Collection

class Parameters(Collection):
    def __init__(self,N=0,prefix=None,names=None,
                 values=None,flags=True,mins=None,maxs=None,
                 weights=None,group='default'):
        self.initialize() # call initialization of a Collection
        #self.__index__ = {}
        self.generate(N=N,prefix=prefix,group=group,
                      names=names,values=values,flags=flags,
                      mins=mins,maxs=maxs,weights=weights)
        self.order = ['group','name','value','min','max','weight','flag','__id__']
                                    
    def __iter__(self):
        return iter(self.getitems()) # strange peculiarity of Python 3
        
    def __len__(self):
        return len(self.ids())
                      
    @property
    def __npar__(self):
        return len(self.ids())
        
    #def __getitem__(self,name):
    #    """
    #    Get parameter object by its name.
    #    """
    #    return self.__dicthash__[self.__index__[name]]
    def __getitem__(self,ID):
        """
        Get parameter object by its name.
        ID = name
        """
        return self.__dicthash__[ID]
        
    def append(self,pars,group='default'):
        """
        Append a parameter to the set. Pars can be either a parameter,
        or a list of parameters.
        """
        if type(pars) is Par:
            pars = [pars]
        elif type(pars) not in {list,tuple}:
            raise Exception('Either Par, or list/tuple of Pars are allowed')
        for par in pars:
            par['group'] = group
            name = par['name']
            #ID = (group,name)
            ID = name
            if ID in self.__dicthash__:
                raise Exception('%s already in __dicthash__ '
                    '(consider renaming parameter)'%str(ID))
            self.__dicthash__[ID] = par
                                       
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
        return np.array(self.getcol('group',IDs=ii))
                    
    def get_names(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array(self.getcol('name',IDs=ii))
        
    def get_values(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array(self.getcol('value',IDs=ii))

    def set_values(self,values,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)        
        for i,val in zip(ii,values):
            self.__dicthash__[i]['value'] = val
        
    def get_bounds(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        data = self.getcols(['min','max'],IDs=ii)
        return list(zip(*data))
        
    def get_weights(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array(self.getcol('weight',IDs=ii))
    
    def get_flags(self,active_only=False,group=None):
        ii = self.get_index(active_only=active_only,group=group)
        return np.array(self.getcol('flag',IDs=ii))
               
    def get_index(self,active_only=False,group=None):   
        cond_group = (lambda v: v['group']==group) if group else (lambda v: True)
        cond_active = (lambda v: v['flag']) if active_only else (lambda v: True)
        cond = lambda v: cond_group(v) and cond_active(v)
        return self.ids(cond)
        
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
                
class Par(dict):
    """
    Fitting parameter storing information about initial values, bounds, rubber weight, and active flag.
    Each parameter has its unique name.
    """
    
    def __init__(self,name,value=None,flag=None,min=None,max=None,weight=None,group=None):
        self['name'] = name
        self['group'] = 'default' if group is None else group
        self['value'] = np.random.random() if value is None else value
        self['flag'] = False if flag is None else flag
        self['min'] = -np.inf if min is None else min
        self['max'] = np.inf if max is None else max
        self['weight'] = 0.0 if weight is None else weight
    
    @property
    def __flag__(self):
        return self['flag']
    
    @property
    def __name__(self):
        return self['name']
    
    def get_value(self):
        return self['value']
        
    def set_value(self,value):
        self['value'] = value
        
    def get_flag(self):
        return self['flag']
        
    def set_flag(self,flag):
        self['flag'] = flag        
        
    def get_bounds(self):
        return [self['min'],self['max']]
        
    def set_bounds(self,bnds):
        self['min'] = bnds[0]
        self['max'] = bnds[1]
        
    def get_min(self):
        return self['min']
        
    def set_min(self,min):
        self['min'] = min

    def get_max(self):
        return self['max']
        
    def set_max(self,max):
        self['max'] = max
        
    def get_weight(self):
        return self['weight']
        
    def set_weight(self,wht):
        self['weight'] = wht
        
    def get_name(self):
        return self['name']    
        
    def get_group(self):
        return self['group']
        
    def __repr__(self):
        group = self.get_group()
        name = self.get_name()
        flag = self.get_flag()
        value = self.get_value()
        min = self.get_min()
        max = self.get_max()
        weight = self.get_weight()
        return '{}:{}:{}({}/[{}:{}])*{}'.format(
            group,name,flag,value,min,max,weight)
