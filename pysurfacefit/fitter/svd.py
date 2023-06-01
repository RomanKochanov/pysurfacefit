import copy

import numpy as np
import numba as nb

from time import time

@nb.njit(cache=True)
def _coeff_mat_numba(X,w,model_func,params,active_indexes):
    """    
    Exploit general PySurfaceFit Model class
    X -> matrix NxM, where N - number of datapoints, M - number of "features"
    w -> weights of datapoints (must have dimensions of Nx1)
    model_func -> function descriptor
    params -> array of model parameters
    active_indexes -> positions of the "active" parameters in this array
    """

    n_points = X.shape[0]
    n_active_pars = len(active_indexes)
        
    mat_ = np.empty((n_points,n_active_pars))
    
    for j in range(n_active_pars): # loop through parameters
        # copy the initial parameters
        params[active_indexes] = np.zeros(n_active_pars)
        # set current parameter to one
        params[active_indexes[j]] = 1
        for i in range(n_points): # loop through datapoints   
            inputs = X[i,:]         
            val = model_func(inputs,params)
            mat_[i,j] = val * w[i,0]
           
    # create offsets to eliminate in the LSQT procedure    
    params[active_indexes] = np.zeros(n_active_pars)
    offset_vect = np.empty((n_points, 1))
    for n in range(n_points):
        inputs = X[n,:]
        offset = model_func(inputs,params) * w[i,0]
        mat_[n,:] -= offset
        offset_vect[n] = offset

    # evaluation of the norm of each column by means of a loop
    scale_vect = np.empty((n_active_pars, 1))
    for n in range(n_active_pars):
        
        # evaluation of the column's norm (stored for later processing)
        col_norm = np.linalg.norm(mat_[:, n])
        scale_vect[n] = col_norm
        
        # scaling the column to unit-length
        mat_[:, n] /= col_norm

    return mat_, scale_vect, offset_vect
    
def _coeff_mat(X,y,w,model):
        
    #if '__initialized_func__' not in model.__dict__ or model.__initialized_func__ is False:
    #    print('=== INITIALIZING FUNC ===')
    #    model.__sympy_initialize_func__()
    #    model.__initialized_func__ = True    
    #    print('=== DONE ===')
        
    model_func = model.__numbified_func__
    
    flags = model.__params__.get_flags(active_only=False)
    active_indexes = np.where(flags)[0]
    
    params = model.__params__.get_values(active_only=False)
    
    print('starting _coeff_mat_numba')
    t = time()
    mat_, scale_vect, offset_vect = _coeff_mat_numba(X,w,model_func,params,active_indexes)
    print('%s sec. elapsed for _coeff_mat_numba'%(time()-t))
    
    y_ = y*w
    
    return mat_, y_, scale_vect, offset_vect
        
@nb.njit(cache=True)
def _fit_x(a, b, scales, offsets):
    # linalg solves ax = b
    det_ = np.linalg.lstsq(a, b+offsets, rcond=-1)[0]
    # due to the stabilization, the coefficients have the wrong scale, which is corrected now
    det_ /= scales

    return det_
