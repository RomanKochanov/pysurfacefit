import copy

from time import time
from scipy import optimize as opt

import numpy as np
from numpy.linalg import pinv, inv

import jeanny3 as j

# default encoder to save numpy types to JSON file
def JSONDefault(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError ("Type %s not serializable" % type(obj))

# ===============================================================================
# ================================== FITTER =====================================
# ===============================================================================
        
class Fitter:
    """
    Class encapsulating the fitting process.
    """

    def __init__(self,model,fitgroups,wpmul=1.0,weighted=True,rubber_on=True,
        method='lm',jac=False,silent=False,**argv):
        """
        model - model function
        inputgrid - grid of inputs
        outputs - function outputs having the same form as grid
        """
        # BASIC PARAMETERS        
        self.__model__ = model # depends on (p,x)
        self.model_final = self.__model__ # MAKE A TEMPORARY ALIAS WHIEN FIRST CREATING THE FITTER OBJECT
        self.__fitgroups__ = fitgroups # groups of points to fit
        self.__method__ = method
        self.__jac__ = jac
        self.__wpmul__ = wpmul
        self.__weighted_fit__ = weighted
        self.__rubber_on__ = rubber_on
        self.__fit_history__ = j.Collection()
        # FIT OPTIONS
        self.__options__ = {} # Nelder-Mead Powell CG BFGS Newton-CG L-BFGS-B TNC COBYLA SLSQP trust-constr dogleg trust-ncg trust-exact trust-krylov
        self.__options__.update(argv)
        self.__icalcfun__ = 0 # function calculation counter
        self.__icalcjac__ = 0 # jacobian calculation counter
        self.__silent__ = silent
        
    def set_options(self,**argv):
        self.__options__ = argv

    def sse(self,p): # simple non-linear general approach to make a target function
        resids = self.residuals(p)
        sum_of_squares = np.sum(resids**2)
        return sum_of_squares
                
    def residuals(self,p): # Levenberg-Marquardt method implemented in the Scipy module
        self.__icalcfun__ += 1
        p = np.array(p); n_p = len(p) 
        self.__model__.__params__.set_values(p,active_only=True) # new format, active=True
        
        # get total number of residuals and create resulting 1D array
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        if self.__rubber_on__ is True:
            n_tot += n_p
            
        resids = np.zeros(n_tot)
        
        # global residual index offset
        offset = 0
               
        # loop through the fit groups
        for grp in self.__fitgroups__:
            if not grp.__active__: continue
            model_calc = self.__model__.calculate(grp.__inputgrid__)
            grp.set_calc_vals(model_calc) # set calculated values for group
            grp.calculate_residuals()
            if self.__weighted_fit__==True:
                resids[offset:(offset+grp.length)] = grp.__weighted_resids__.flatten()
            else:
                resids[offset:(offset+grp.length)] = grp.__unweighted_resids__.flatten()
            offset += grp.length
               
        # do "rubber"
        p_resids = p-self.__params_initial__.get_values(active_only=True) # unweighted
        resids_param = self.__wpmul__*self.__model__.__wp__*p_resids # weighted
        if self.__rubber_on__ is True:                        
            resids[offset:(offset+n_p)] = resids_param
            offset += n_p
        
        # output iteration statistics        
        sum_of_squares_tot = np.sum(resids**2)
        # sum_of_squares_grp = np.sum(resids[n_p:]**2); print('sum_of_squares_grp>>>','%.9e'%sum_of_squares_grp)
        sum_of_squares_rub = np.sum(resids_param**2)
        if not self.__silent__ and self.__icalcfun__%1==0:
            print('CALC FUN %5d>>>  DIST:%13.9e  SSE_RUB:%13.9e SSE_TOT:%13.9e %20s'%\
                 (self.__icalcfun__,np.sqrt(np.sum(p_resids**2)),
                  sum_of_squares_rub,sum_of_squares_tot,
                  '==> WEIGHTED_FIT' if self.__weighted_fit__ else '==> UNWEINGHTED_FIT'
                  ))
        if not self.__silent__: print('=====data group statistics=====')
        stat = self.__fitgroups__.get_group_stats()
        if not self.__silent__: stat.tabulate(floatfmt=[
            None,   # GROUP 
            '5.0f', # N   
            '5.1e', # MIN_WHT 
            '5.1e', # MAX_WHT  
            '7.3e', # MIN_WHT_RES, 
            '7.3e', # MAX_WHT_RES,  
            '7.3e', # MIN_UNWHT_RES, 
            '7.3e', # MAX_UNWHT_RES
            '7.3e', # WHT_SSE  
            '7.3e', # UNWHT_SSE 
            '7.3e', # WHT_SD   
            '7.3e', # UNWHT_SD 
        ])
        if not self.__silent__: print('===============================')
                
        # Save history item to fit history collection.
        history_item = {}
        history_item['stat_dicthash'] = copy.deepcopy(stat.__dicthash__)
        history_item['params_dicthash'] = copy.deepcopy(self.__model__.__params__.__dicthash__)
        self.__fit_history__.update(history_item)
                
        return resids
        
    def residuals_jac(self,p):
        self.__icalcjac__ += 1
        p = np.array(p); n_p = len(p)        
        self.__model__.__params__.set_values(p,active_only=True) # new format, active=True
        
        # get total number of residuals and create resulting 2D matrix
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        if self.__rubber_on__ is True:
            n_tot += n_p
            
        resids = np.zeros([n_tot,n_p])
        
        # global residual index offset
        offset = 0
                
        # loop through the fit groups
        for grp in self.__fitgroups__:
            if not grp.__active__: continue
            model_calc_jac = self.__model__.calculate_jac(grp.__inputgrid__) # calculate Jacobian, flatten the initial grid structure
            model_calc = self.__model__.calculate(grp.__inputgrid__) # this step can be omitted to gain some (presumably low) speed-up
            grp.set_calc_vals(model_calc) # set calculated values for group; can be optimized as well by "caching"
            grp.calculate_residual_derivs()
            if self.__weighted_fit__==True:
                resids[offset:(offset+grp.length)] = model_calc_jac * np.reshape(grp.__weighted_resid_derivs__.flatten(),[grp.length,1])
            else:
                resids[offset:(offset+grp.length)] = model_calc_jac * np.reshape(grp.__unweighted_resid_derivs__.flatten(),[grp.length,1])
            offset += grp.length        

        # do "rubber"        
        if self.__rubber_on__ is True:                  
            resids_param = self.__wpmul__*self.__model__.__wp__*np.eye(n_p) # weighted derivatives
            resids[offset:(offset+n_p)] = resids_param
            offset += n_p       

        # output iteration statistics
        if self.__icalcjac__%1==0:
            p_resids = p-self.__params_initial__.get_values(active_only=True)
            if not self.__silent__: print('CALC JAC %5d>>>  DIST:%13.9e'%(self.__icalcjac__,np.sqrt(np.sum(p_resids**2))))

        return resids
        
    def set_weights(self,wmod):
        self.__wmod__ = np.array(wmod)
        
    #def set_penalty(self,penalty_inputgrid,penalty_output,wpenmul):
    #    self.__penalty_inputgrid__ = penalty_inputgrid
    #    self.__penalty_output__ = np.array(penalty_output)
    #    self.__wpen__ = np.ones(penalty_inputgrid.__shape__)
    #    self.__wpenmul__ = wpenmul

    def fit_minimize(self): # a family of methods from scipy.optimize.minimize
        if not self.__silent__: print('USING SCIPY.OPTIMIZE.MINIMIZE: METHOD=',self.__method__)
        argv = self.__options__.copy()
        argv['bounds'] = self.__model__.__params__.get_bounds(active_only=True)
        if self.__jac__:
            #argv['jac'] = lambda p:self.residuals_jac(p)
            raise NotImplementedError # need to implement the derivatives for minimize first
        res = opt.minimize(self.sse,self.__model__.__p__,method=self.__method__,**argv)
        self.__fitgroups__.split_global_jacobian(res.jac) # split Jacobian into groups to make statistics
        return res
    
    def fit_basinhopping(self):
        if not self.__silent__: print('USING SCIPY.OPTIMIZE.BASINHOPPING/ANNEAL: METHOD=',self.__method__)
        if not self.__silent__: print('WARNING: PARAMETER BOUNDS ARE IGNORED')
        argv = self.__options__.copy()
        del argv['max_nfev']
        bounds = self.__model__.__params__.get_bounds(active_only=True)
        bounds = list(zip(*bounds))
        class MyBounds(object): #https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html
            def __init__(self, xmax=bounds[0], xmin=bounds[1]):
                self.xmax = np.array(xmax)
                self.xmin = np.array(xmin)
            def __call__(self, **kwargs):
                x = kwargs["x_new"]
                tmax = bool(np.all(x <= self.xmax))
                tmin = bool(np.all(x >= self.xmin))
                return tmax and tmin     
        argv['accept_test'] = MyBounds()
        #if self.__jac__:
        #    #argv['jac'] = lambda p:self.residuals_jac(p)
        #    raise NotImplementedError # need to implement the derivatives for minimize first
        res = opt.basinhopping(self.sse,self.__model__.__p__,**argv)
        return res

    def fit_least_squares(self): # Least squares method
        if not self.__silent__: print('USING SCIPY.OPTIMIZE.LEAST_SQUARES: METHOD=',self.__method__)
        argv = self.__options__.copy()
        if self.__method__ in {'trf','dogbox'}:
            argv['bounds'] = list(zip(*self.__model__.__params__.get_bounds(active_only=True)))
        elif self.__method__ == 'lm':
            if not self.__silent__: print('Warning: skipping bounds for the "lm" method')
        if self.__jac__:
            argv['jac'] = lambda p:self.residuals_jac(p)
        if not self.__silent__: print('METHOD OPTIONS:',argv)
        res = opt.least_squares(self.residuals,self.__model__.__p__,
            method=self.__method__,**argv) # with user-supplied jacobian
        #self.__fitgroups__.split_global_jacobian(res.jac) # split Jacobian into groups to make statistics
        #self.calculate_covariance_matrix(res) # calculate covariance matrix
        # save Jacobian on disc
        jac_file = self.__model__.__class__.__name__+'.jac'
        if not self.__silent__: 
            print('Saving Jacobian matrix to %s.npy'%jac_file)
            np.save(jac_file,res.jac)
        return res
        
    def fit(self):        
        t = time()
        self.__params_initial__ = copy.deepcopy(self.__model__.__params__)
        icalcfun_old = self.__icalcfun__
        icalcjac_old = self.__icalcjac__
        if not self.__silent__: print('BEGIN FIT')
        if self.__method__ in {'lm','trf','dogbox'}:
            res = self.fit_least_squares()
        elif self.__method__ in {'Nelder-Mead','Powell','CG','BFGS',
                                 'Newton-CG','L-BFGS-B','TNC','COBYLA',
                                 'SLSQP','trust-constr','dogleg',
                                 'trust-ncg','trust-exact','trust-krylov'}:
            res = self.fit_minimize()
        elif self.__method__ in {'basinhopping','anneal'}:
            res = self.fit_basinhopping()
        else:
            raise Exception('unknown method: ',self.__method__)
        #self.__params_final__ = copy.deepcopy(self.__model__.__params__)        
        if not self.__silent__: print('END FIT')
        if not self.__silent__: print('%f seconds elapsed, %d func evals, %d jac evals'%\
            (time()-t,self.__icalcfun__-icalcfun_old,self.__icalcjac__-icalcjac_old))
        self.__result__ = res
        # create model final
        self.model_final = self.__model__ # create a link instead of a separate object
        # save fit history
        if not self.__silent__: 
            self.__fit_history__.export_json_list(
                self.__model__.__class__.__name__+'.fit_history',default=JSONDefault)
        
    @property
    def model_initial(self):
        if '__dll_func__' in self.__model__.__dict__:
            self.__model__.__dll_func__ = None # make deepcopy work
        cpy = copy.deepcopy(self.__model__)
        cpy.__params__ = self.__params_initial__
        return cpy

    #@property
    #def model_final(self,*inputs):
    #    if '__dll_func__' in self.__model__.__dict__:
    #        self.__model__.__dll_func__ = None # make deepcopy work
    #    cpy = copy.deepcopy(self.__model__)
    #    #cpy.__params__ = self.__params_final__        
    #    return cpy
        
    #def calculate_covariance_matrix(self,fit_result):
    #    """Calculate covariance CovB between parameters.
    #       https://stats.stackexchange.com/questions/231868/relation-between-covariance-matrix-and-jacobian-in-nonlinear-least-squares"""
    #    #n = len(fit_result.fun)
    #    #X = fit_result.jac
    #    n = len(fit_result.fun)
    #    X = fit_result.jac
    #    MSE = np.sum(fit_result.fun**2)/n
    #    CovB = pinv(X.T @ X)*MSE
    #    #CovB = inv(X.T @ X)*MSE
    #    self.__covb__ = CovB
                
    def __repr__(self):
        if '__result__' in self.__dict__:
            return str(self.__result__)
        else:
            return '<no __result__>'
