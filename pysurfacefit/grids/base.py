from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import pylab as pl
import numpy as np

def assert_(desc,vars):
    assert set(vars) not in set(desc.keys()),\
            'name mismatch'

def flatten(*args):
    if len(args)==1:
        return args[0].flatten()
    else:
        return [_.flatten() for _ in args]
        
class Grid:
    """
    desc: {'v1':grid1, 'v2':grid2 ...}
    vars: ['v1','v2',...]
    """
    def __init__(self,*argc):
        desc = dict([(argc[i],argc[i+1]) for i in range(0,len(argc),2)]); vars = argc[0::2]
        assert_(desc,vars)
        self.__desc__ = {_:np.array(desc[_]) for _ in desc}
        self.__vars__ = vars
        grids = [self.__desc__[name] for name in vars]        
        meshes = np.meshgrid(*grids)
        self.__meshes__ = meshes
        self.__size__ = meshes[0].size
        self.__shape__ = meshes[0].shape
        self.__hash__ = {name:i for i,name in enumerate(vars)} # hash for lookup mesh by the variable name
    
    def __getitem__(self,varname):
        return self.__meshes__[self.__hash__[varname]]
        
    def flatten(self):
        return flatten(*self.__meshes__)
        
    def get_meshes(self,ind=None,flat=False,squeeze=True):
        # https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.squeeze.html
        meshes = self.__meshes__
        if flat:
            meshes = flatten(*meshes)
        if squeeze:
            meshes = [np.squeeze(_) for _ in meshes]
        if len(meshes)==1: meshes=meshes[0]
        if ind is None:
            return meshes
        else:
            return [mesh[ind] for mesh in meshes]
                                    
    def calculate(self,func,flat=False,squeeze=True,dtype=np.float64): # calculate arbitrary function on the grid; SEEMS IT CAN BE NUMBIFIED! (ndenumerate IS SUPPORTED BY NUMBA)        
        meshes = self.__meshes__
        res = np.empty(meshes[0].shape,dtype=dtype)
        for i,_ in np.ndenumerate(res): # np.ndenumerate, np.nditer
            arg = [mesh[i] for mesh in meshes]
            res[i] = func(*arg)
        if flat:
            res = flatten(res)
        if squeeze:
            res = np.squeeze(res)
        return res
        
class List:
    """
    desc: {'v1':grid1, 'v2':grid2 ...}
    vars: ['v1','v2',...]
    """
    def __init__(self,*argc):
        desc = dict([(argc[i],argc[i+1]) for i in range(0,len(argc),2)]); vars = argc[0::2]
        assert_(desc,vars)
        self.__desc__ = {_:np.array(desc[_]) for _ in desc}
        self.__vars__ = vars
        meshes = [self.__desc__[name] for name in vars]
        self.__meshes__ = meshes
        self.__size__ = meshes[0].size
        self.__shape__ = meshes[0].shape
        self.__hash__ = {name:i for i,name in enumerate(vars)} # hash for lookup mesh by the variable name
        
    def __getitem__(self,varname):
        return self.__meshes__[self.__hash__[varname]]        
        
    def flatten(self): # newly added
        return flatten(*self.__meshes__)
                        
    def get_meshes(self,ind=None,flat=False):
        meshes = self.__meshes__
        if len(meshes)==1: meshes=meshes[0]
        if ind is None:
            return meshes
        else:
            return [mesh[ind] for mesh in meshes]
                                    
    def calculate(self,func,flat=False,dtype=np.float64): # calculate arbitrary function on the grid; SEEMS IT CAN BE NUMBIFIED! (ndenumerate IS SUPPORTED BY NUMBA)
        meshes = self.__meshes__
        res = np.empty(meshes[0].shape,dtype=dtype)
        for i,_ in np.ndenumerate(res): # np.ndenumerate, np.nditer
            arg = [mesh[i] for mesh in meshes]            
            res[i] = func(*arg)
        return res