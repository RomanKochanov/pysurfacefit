import numpy as np

from . import fortran as modfort 

from .base import Model

from ..fitpars import Par, Parameters

#from numba import prange,njit,jitclass,int32,float32 # jitclass doesn't work with abstract methods
from numba import prange,njit,int32,float32 # jitclass doesn't work with abstract methods

# IDEA: make a Parameter class, which supports bounds etc...
# this will help to get rid of "__npar__"
# Drawback: need to re-define all basic functions and operators
# for Parameter class
#https://www.python-course.eu/python3_magic_methods.php
         
class PolynomialModel1D(Model):
    """
    1D polynomial model 
    """
    def __init__(self,maxpower):
        self.__params__ = Parameters(maxpower)

    def __func__(self,params,x):
        p = params.get_values()
        return np.polyval(p,x)
            
def poly2d(p,x,y,powers):
    n = len(p)
    sum = 0.0
    for i in range(n):
        sum += p[i]*x**powers[i][0]*y**powers[i][1]
    return sum
fast_poly2d = njit(poly2d)

def generate_powers(maxpower):
    powers = []
    for pow_x in range(maxpower+1):
        for pow_y in range(maxpower+1-pow_x):
            powers.append([pow_x,pow_y])
    return powers

def generate_powers_layer_3d(total_power):
    """
    Generate layer of powers for 3d polynomial for a given total power (unsymmetrized)
    """
    powers = []
    for j_plus_k in range(total_power,-1,-1):
        i = total_power - j_plus_k
        for j in range(j_plus_k+1):
            k = j_plus_k - j
            powers.append([i,j,k])
    return powers

def generate_powers_3d_rsymm(maxpower):
    """
    Generate powers for polynomial symmetrized by first 2 coordinates (r1,r2)
    """
    powers = []
    for pow_alpha in range(maxpower+1): # power of alpha
        for pow_r1_plus_r2 in range(maxpower-pow_alpha+1): # sum of r1 and r2 powers
            for pow_r1 in range(pow_r1_plus_r2+1): # power of r1
                pow_r2 = pow_r1_plus_r2-pow_r1
                if pow_r1<pow_r2: continue
                powers.append([pow_r1,pow_r2,pow_alpha])
    return powers

def poly3d(p,x,y,z,powers):
    n = len(p)
    sum = 0.0
    for i in range(n):
        sum += p[i]*x**powers[i][0]*y**powers[i][1]*z**powers[i][2]
    return sum
fast_poly3d = njit(poly2d)
    
def poly3d_rsymm(p,r1,r2,alpha,powers):
    """
    3D polynomial symmetrized by first 2 coordinates (r1,r2)
    """
    n = len(p)
    sum = 0.0
    for i in range(n):
        pow_r1 = powers[i][0]
        pow_r2 = powers[i][1]
        pow_alpha = powers[i][2]
        sum += p[i]*(r1**pow_r1*r2**pow_r2+r1**pow_r2*r2**pow_r1)*alpha**pow_alpha
    return sum
fast_poly3d_rsymm = njit(poly3d_rsymm)
    
class PolynomialModel2D(Model):
    """
    2D polynomial model
    """
    def __init__(self,maxpower,fast=False):
        self.__powers__ = np.array(generate_powers(maxpower))
        self.__params__ = Parameters(len(self.__powers__),
            names=['%d_%d'%tuple(_) for _ in self.__powers__],
            group='linear')
        self.__fast__ = fast
                    
    def __func__(self,params,x,y):
        p = params.get_values()
        if self.__fast__:
            res = fast_poly2d(p,x,y,self.__powers__)
        else:
            res = poly2d(p,x,y,self.__powers__)
        return res

# ========================================        
# ========================================        
# MODELS FOR FITTING THE 666 ISOTOPOLOGUE        
# ========================================        
# ========================================        
        
class DBOCModel_1(Model):
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            [0,     0], # 0
            
            [1,     0], # 1
            [0,     1], # 1
            
            [1,     1], # 2
            [0,     2], # 2
            [2,     0], # 2
            
            [3,     0], # 3
            [2,     1], # 3
            [1,     2], # 3
            [0,     3], # 3
            
            [4,     0], # 4
            [3,     1], # 4
            [2,     2], # 4
            [1,     3], # 4
            [0,     4], # 4            
        ])
        self.__params__ = Parameters(len(self.__powers__),
            names=['%d_%d'%tuple(_) for _ in self.__powers__])
            
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        p = params.get_values()
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
        if self.__fast__:
            res = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            res = poly2d(p,r_part,alpha_part,self.__powers__)
        return res
        
class DBOCModel_2(Model): # GOOD MODEL TO FIT DBOC(6,6) FOR 666 ISOTOPOLOGUE; DON'T CHANGE
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
        #   r_symm  alpha
            [0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
        ])
                
        # polynomial parameters, active by default, no bounds
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,flag=True)
        
        # gaussian parameters, with bounds
        self.__params__.append(group='gaussian', pars=[
            Par(name='r_gauss', value=2.01, min=2.0, max=4.0,  flag=True),
            Par(name='a_gauss', value=9.99964, min=2.0, max=10.0, flag=True),
        ])
        
        # structure parameters, with bounds
        weight = 0.0
        self.__params__.append(group='structure', pars=[
            Par(name='r0_struct',value=2.95988, min=1.0, max=5.0,  weight=weight, flag=True),
            Par(name='a_struct', value=2.1391,  min=1,   max=50.0, weight=weight, flag=True),
            Par(name='h_struct', value=7.0, min=-np.inf, max=np.inf, weight=weight, flag=True),
            Par(name='apoly0_struct', value=3.01235,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly1_struct', value=1.06193,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly2_struct', value=0.269482, min=-np.inf, weight=weight, max=np.inf, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=0.0,  min=-30.0,   max=30.0,  flag=False),
        ])
                
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        r_gauss = params['r_gauss'].get_value()
        a_gauss = params['a_gauss'].get_value()
        gauss_r1 = np.exp(-a_gauss*(r1_au-re_au)**2)
        gauss_r2 = np.exp(-a_gauss*(r2_au-re_au)**2)
        gauss = gauss_r1*gauss_r2
        
        # get structural parameters
        r0_struct = params['r0_struct'].get_value()
        a_struct = params['a_struct'].get_value()
        h_struct = params['h_struct'].get_value()
        apoly0_struct = params['apoly0_struct'].get_value()
        apoly1_struct = params['apoly1_struct'].get_value()
        apoly2_struct = params['apoly2_struct'].get_value()
        struct_r1 = 0.5*(np.tanh(a_struct*(r1_au-r0_struct)/r0_struct)+1)
        struct_r2 = 0.5*(np.tanh(a_struct*(r2_au-r0_struct)/r0_struct)+1)
        struct = h_struct*struct_r1*struct_r2
        struct *= (apoly0_struct+apoly1_struct*alpha_part+apoly2_struct*alpha_part**2)
        
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # get final result
        res = poly*gauss+struct+global_shift
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
        
        return res
      
class DBOCModel_test(Model): # SANDBOX MODEL FOR TESTING PURPOSES
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            [0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
            
            [3,     0], # 3
            [2,     1], # 3
            [1,     2], # 3
            [0,     3], # 3
            
            [4,     0], # 4
            [3,     1], # 4
            [2,     2], # 4
            [1,     3], # 4
            [0,     4], # 4            
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 5.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='r_gauss', value=2.5, min=2.0, max=3.0, weight=weight, flag=True),
            Par(name='a_gauss', value=3.5, min=2.0, max=10.0, weight=weight, flag=True),
        ])
        
        # structure parameters, with bounds
        weight = 0.0
        self.__params__.append(group='structure', pars=[
            Par(name='r0_struct',value=2.95988, min=1.0, max=20.0,  weight=weight, flag=True),
            Par(name='a_struct', value=2.1391,  min=1,   max=50.0, weight=weight, flag=True),
            Par(name='h_struct', value=7.0, min=-np.inf, max=np.inf, weight=weight, flag=True),
            Par(name='r0_apoly1_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='r0_apoly2_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='r0_apoly3_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='a_apoly1_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='a_apoly2_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='a_apoly3_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly0_struct', value=1.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly1_struct', value=0.0,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly2_struct', value=0.0, min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly3_struct', value=0.0, min=-np.inf, weight=weight, max=np.inf, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=0.0,  min=-30.0,   max=30.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        r_gauss = params['r_gauss'].get_value()
        a_gauss = params['a_gauss'].get_value()
        gauss_r1 = np.exp(-a_gauss*(r1_au-r_gauss)**2)
        gauss_r2 = np.exp(-a_gauss*(r2_au-r_gauss)**2)
        gauss = gauss_r1*gauss_r2
                
        # get structural parameters
        r0_struct = params['r0_struct'].get_value()
        a_struct = params['a_struct'].get_value()
        h_struct = params['h_struct'].get_value()
        r0_apoly1_struct = params['r0_apoly1_struct'].get_value()
        r0_apoly2_struct = params['r0_apoly2_struct'].get_value()
        r0_apoly3_struct = params['r0_apoly3_struct'].get_value()
        a_apoly1_struct = params['a_apoly1_struct'].get_value()
        a_apoly2_struct = params['a_apoly2_struct'].get_value()
        a_apoly3_struct = params['a_apoly3_struct'].get_value()
        apoly0_struct = params['apoly0_struct'].get_value()
        apoly1_struct = params['apoly1_struct'].get_value()
        apoly2_struct = params['apoly2_struct'].get_value()
        apoly3_struct = params['apoly3_struct'].get_value()
        r0_struct = r0_struct*(1+r0_apoly1_struct*alpha_part+\
                                 r0_apoly2_struct*alpha_part**2+\
                                 r0_apoly3_struct*alpha_part**3)
        a_struct = a_struct*(1+a_apoly1_struct*alpha_part+\
                               a_apoly2_struct*alpha_part**2+\
                               a_apoly3_struct*alpha_part**3)
        #struct_r1 = 1/2*(np.tanh(a_struct*(r1_au-r0_struct)/r0_struct)+1)
        #struct_r2 = 1/2*(np.tanh(a_struct*(r2_au-r0_struct)/r0_struct)+1)
        struct_r1 = 1/2*(np.tanh(a_struct*(r1_au-r0_struct))+1)
        struct_r2 = 1/2*(np.tanh(a_struct*(r2_au-r0_struct))+1)
        struct = h_struct*struct_r1*struct_r2
        struct *= (apoly0_struct+apoly1_struct*alpha_part+\
                   apoly2_struct*alpha_part**2+apoly3_struct*alpha_part**3)
                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # get final result
        res = poly*gauss+struct+global_shift
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res

class DBOCModel_666_22(Model): # MODEL TO FIT DBOC(2,2) FOR 666 ISOTOPOLOGUE
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            #[0,     2], # 2   
            [2,     0], # 2   
            
            [3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=1.5, min=-50.0, max=50.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=0.5, min=0.4, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=107.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=35.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=10.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res

class DBOCModel_666_44(Model): # MODEL TO FIT DBOC(4,4) FOR 666 ISOTOPOLOGUE
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            #[0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=1.5, min=-50.0, max=50.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=0.5, min=0.4, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=107.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=35.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=10.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res

class DBOCModel_666_66(Model): # LAME MODEL TO FIT DBOC(6,6) FOR 666 ISOTOPOLOGUE
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=2.5, min=2.0, max=3.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=10.0, min=5.0, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=107.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=20.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=10.0,  min=5.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res

class DBOCModel_666_66_2(Model): # model to refit the full DBOC correction for 66; direct copy of DBOCModel_666_66 with some initial parameter changes
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=2.4, min=2.0, max=3.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=10.0, min=5.0, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=117.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=15.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=5.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res
        
class DBOCModel_test2(Model): # SANDBOX MODEL FOR TESTING PURPOSES
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            #[0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=1.5, min=-50.0, max=50.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=0.5, min=0.4, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=107.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=35.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=10.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = poly*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res
        
class DBOCModel_3(Model):
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
        #   r_symm  alpha
            [0,     0], # 0
            
            [1,     0], # 1
            [0,     1], # 1
            
            [1,     1], # 2
            [0,     2], # 2
            [2,     0], # 2
        ])
                
        # polynomial parameters, active by default, no bounds
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,flag=True)
        
        # gaussian parameters, with bounds
        self.__params__.append(group='gaussian', pars=[
            Par(name='r_gauss', value=2.01, min=2.0, max=4.0,  flag=True),
            Par(name='a_gauss', value=9.99964, min=2.0, max=10.0, flag=True),
        ])
        
        # structure parameters, with bounds
        weight = 10.0
        self.__params__.append(group='structure', pars=[
            Par(name='r0_struct',value=2.95988, min=1.0, max=5.0,  weight=weight, flag=True),
            Par(name='a_struct', value=2.1391,  min=1,   max=50.0, weight=weight, flag=True),
            Par(name='h_struct', value=7.0, min=-np.inf, max=np.inf, weight=weight, flag=False),
            Par(name='apoly0_struct', value=3.01235,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly1_struct', value=1.06193,  min=-np.inf, weight=weight, max=np.inf, flag=True),
            Par(name='apoly2_struct', value=0.269482, min=-np.inf, weight=weight, max=np.inf, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=0.0,  min=-30.0,   max=30.0,  flag=False),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        r_gauss = params['r_gauss'].get_value()
        a_gauss = params['a_gauss'].get_value()
        gauss_r1 = np.exp(-a_gauss*(r1_au-re_au)**2)
        gauss_r2 = np.exp(-a_gauss*(r2_au-re_au)**2)
        gauss = gauss_r1*gauss_r2
        
        # get structural parameters
        r0_struct = params['r0_struct'].get_value()
        a_struct = params['a_struct'].get_value()
        h_struct = params['h_struct'].get_value()
        apoly0_struct = params['apoly0_struct'].get_value()
        apoly1_struct = params['apoly1_struct'].get_value()
        apoly2_struct = params['apoly2_struct'].get_value()
        struct_r1 = 1/2*(np.tanh(a_struct*(r1_au-r0_struct)/r0_struct)+1)
        struct_r2 = 1/2*(np.tanh(a_struct*(r2_au-r0_struct)/r0_struct)+1)
        struct = h_struct*struct_r1*struct_r2
        struct *= (apoly0_struct+apoly1_struct*alpha_part+apoly2_struct*alpha_part**2)
        
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # get final result
        res = poly*gauss+struct+global_shift
        
        return res


        
class TestModelFortran(modfort.ModelFortran):
    """
    Sample Fortran model.
    """
    def __init__(self):
        self.__name__ = 'my_test_1'
        self.__maxpower__ = 10
        self.__powers__ = np.array(generate_powers(self.__maxpower__))
        self.__params__ = Parameters(N=len(self.__powers__))
        self.__data__ = modfort.Data(
            modfort.FortranValue('Shift','double precision',100.0),
            modfort.FortranMatrix('Powers','integer',self.__powers__),
            )
        self.__code__ = """
        subroutine foo(res,p,x,y)
        integer :: i
        !res = p(1)*x+p(2)*y+p(3)*x**2+p(4)*y**2+Shift
        res = 0.0d0
        do i=1,np__
            res = res + p(i)*x**Powers(i,1)*y**Powers(i,2)
        end do
        !print *,'p',p
        !print *,'x',x,'y',y,'res',res
        end subroutine
        """
        super().__init__()
            
class PolynomialModelFortran():
    pass
    

# ========================================        
# ========================================        
# MODELS FOR FITTING THE 888 ISOTOPOLOGUE        
# ========================================        
# ========================================        

class DBOCModel_888_66_2(Model): # model to refit the full DBOC correction for 66; direct copy of DBOCModel_666_66_2
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=2.4, min=2.0, max=3.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=10.0, min=5.0, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=117.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=15.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=5.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = (1+poly)*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res
        
# ========================================        
# ========================================        
# MODELS FOR FITTING THE 777 ISOTOPOLOGUE        
# ========================================        
# ========================================        
        
class DBOCModel_777_66_2(Model): # model to refit the full DBOC correction for 66, iso 777; direct copy of DBOCModel_666_66_2
    """
    Model for DBOC points fitting
    """
    def __init__(self,fast=False):
        self.__components__ = {}
        self.__fast__ = fast
        self.__powers__ = np.array([
        #   r_symm  alpha
            #[0,     0], # 0   
            
            [1,     0], # 1   
            [0,     1], # 1   
            
            [1,     1], # 2   
            [0,     2], # 2   
            [2,     0], # 2   
            
            #[3,     0], # 3
            #[2,     1], # 3
            #[1,     2], # 3
            #[0,     3], # 3
            
            #[4,     0], # 4
            #[3,     1], # 4
            #[2,     2], # 4
            #[1,     3], # 4
            #[0,     4], # 4
        ])
                
        # polynomial parameters, active by default, no bounds
        weight = 0.0
        self.__params__ = Parameters(group='polynom',N=len(self.__powers__),
            prefix='lin_',names=['%d_%d'%tuple(_) for _ in self.__powers__],
            value=0.0,weight=weight,flag=True)
        
        # gaussian parameters, with bounds
        weight = 0.0
        self.__params__.append(group='gaussian', pars=[
            Par(name='center_rgauss', value=2.4, min=2.0, max=3.0, weight=weight, flag=True),
            Par(name='slope_rgauss', value=10.0, min=5.0, max=20.0, weight=weight, flag=True),
            Par(name='center_agauss', value=117.0, min=-180.0, max=180.0, weight=weight, flag=True),
            Par(name='slope_agauss', value=4.1, min=2.0, max=10.0, weight=weight, flag=True),
            Par(name='mul_gauss', value=15.0, min=0.0, max=100.0, weight=weight, flag=True),
        ])

        # global shift parameters, with bounds
        self.__params__.append(group='globshift', pars=[
            Par(name='global_shift', value=15.0,  min=5.0,   max=20.0,  flag=True),
        ])
        
    def __func__(self,params,r1_au,r2_au,alpha_deg):
        # common variables
        re_au = 2.4
        alpha_rad = alpha_deg*np.pi/180.
        alphae_rad = 117.0*np.pi/180.
        r_part = (r1_au+r2_au)/re_au
        alpha_part = (np.cos(alpha_rad)-\
                np.cos(alphae_rad))/np.cos(alphae_rad)
                                
        # get polynomial part
        p = params.get_values(group='polynom')
        if self.__fast__:
            poly = fast_poly2d(p,r_part,alpha_part,self.__powers__)
        else:
            poly = poly2d(p,r_part,alpha_part,self.__powers__)
            
        # get gaussian parameters
        center_rgauss = params['center_rgauss'].get_value()
        slope_rgauss = params['slope_rgauss'].get_value()
        center_agauss = params['center_agauss'].get_value()
        slope_agauss = params['slope_agauss'].get_value()
        mul_gauss = params['mul_gauss'].get_value()
        gauss_r1 = np.exp(-slope_rgauss*((r1_au-center_rgauss)/center_rgauss)**2)
        gauss_r2 = np.exp(-slope_rgauss*((r2_au-center_rgauss)/center_rgauss)**2)
        gauss_alpha = np.exp(-slope_agauss*((alpha_deg-center_agauss)/center_agauss)**2)
        #gauss_alpha = 1
        gauss = mul_gauss*gauss_r1*gauss_r2*gauss_alpha
                                           
        # get global shift
        global_shift = params['global_shift'].get_value()
        
        # structure is set to zero
        struct = 0.0
        
        # get final result
        res = global_shift-(1+poly)*gauss
        
        # set components
        self.__components__['gauss'] = gauss
        self.__components__['poly'] = poly
        self.__components__['polygauss'] = (1+poly)*gauss
        self.__components__['struct'] = struct 
        self.__components__['globshift'] = global_shift
            
        return res
