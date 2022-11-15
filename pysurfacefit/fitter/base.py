import copy

from time import time
from scipy import optimize as opt

import numpy as np
from numpy.linalg import pinv, inv

# ===============================================================================
# ================================== FITTER =====================================
# ===============================================================================
        
class Fitter:
#    def __init__(self,model,inputgrid,output,
#            penalty_inputgrid=None,penalty_output=[],
#            penalty_upper_on=False,penalty_lower_on=False,
#            wpen=None,wpenmul=1.0,
#            method='LM',wmod=None,wpmul=1.0,jac=False,**argv):
    def __init__(self,model,fitgroups,wpmul=1.0,weighted=True,rubber_on=True,method='LM',jac=False,**argv):
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
        # FIT OPTIONS
        self.__options__ = {} # Nelder-Mead Powell CG BFGS Newton-CG L-BFGS-B TNC COBYLA SLSQP trust-constr dogleg trust-ncg trust-exact trust-krylov
        self.__options__.update(argv)
        self.__icalcfun__ = 0 # function calculation counter
        self.__icalcjac__ = 0 # jacobian calculation counter
        
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
        
        print('ENTERING RESIDUALS')
        #print('len(p)>>>',len(p))
               
        # get total number of residuals and create resulting 1D array
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        if self.__rubber_on__ is True:
            n_tot += n_p
            
        resids = np.zeros(n_tot)
        
        # global residual index offset
        offset = 0
        
        #print('self.__rubber_on__>>>',self.__rubber_on__)
        #print('self.__rubber_on__ is True>>>',self.__rubber_on__ is True)
        
        # do "rubber" first
        p_resids = p-self.__params_initial__.get_values(active_only=True) # unweighted
        resids_param = self.__wpmul__*self.__model__.__wp__*p_resids # weighted
        if self.__rubber_on__ is True:                        
            resids[offset:(offset+n_p)] = resids_param
            #print('offset,offset+n_p>>>',offset,offset+n_p)
            offset += n_p
            
        print('RESIDUALS: loop through the fit groups')
        
        # loop through the fit groups
        for grp in self.__fitgroups__:
            #print('RESIDUALS: grp.__name__>>>',grp.__name__)
            if not grp.__active__: continue
            #print('RESIDUALS: start calc on grid...')
            model_calc = self.__model__.calculate(grp.__inputgrid__)
            #print('RESIDUALS: ...done')
            # DEBUG
            #m1,m2,m3 = grp.__inputgrid__.get_meshes(flat=True)
            #print(np.array(list(zip(m1,m2,m3,model_calc,grp.__output__))))
            #print('model_calc>>>',model_calc)
            # //DEBUG
            #print('RESIDUALS: setting calc vals and calculating residuals...')
            grp.set_calc_vals(model_calc) # set calculated values for group
            grp.calculate_residuals()
            #print('RESIDUALS: ...done')
            if self.__weighted_fit__==True:
                resids[offset:(offset+grp.length)] = grp.__weighted_resids__.flatten()
            else:
                resids[offset:(offset+grp.length)] = grp.__unweighted_resids__.flatten()
            #print('offset,offset+grp.length>>>',offset,offset+grp.length)
            offset += grp.length
        
        #print('self.__weighted_fit__>>>',self.__weighted_fit__)
        #print('np.max(np.abs(resids))>>>',np.max(np.abs(resids)))
        
        #print('RESIDUALS: print statistics')
        
        # output iteration statistics        
        sum_of_squares_tot = np.sum(resids**2)
        # sum_of_squares_grp = np.sum(resids[n_p:]**2); print('sum_of_squares_grp>>>','%.9e'%sum_of_squares_grp)
        sum_of_squares_rub = np.sum(resids_param**2)
        if self.__icalcfun__%1==0:
            print('CALC FUN %5d>>>  DIST:%13.9e  SSE_RUB:%13.9e SSE_TOT:%13.9e %20s'%\
                 (self.__icalcfun__,np.sqrt(np.sum(p_resids**2)),
                  sum_of_squares_rub,sum_of_squares_tot,
                  '==> WEIGHTED_FIT' if self.__weighted_fit__ else '==> UNWEINGHTED_FIT'
                  ))
        print('=====data group statistics=====')
        stat = self.__fitgroups__.get_group_stats()
        stat.tabulate(floatfmt=[
        #   GROUP  N       MIN_WHT  MAX_WHT  WHT_SSE  UNWHT_SSE  WHT_SD   UNWHT_SD
            None, '5.0f', '5.1e',   '5.1e',  '7.3e',  '7.3e',    '7.3e',  '7.3e'])
        print('===============================')
                
        #print('QUITTING RESIDUALS')
                
        return resids
        
    def residuals_jac(self,p):   # SHOULD BE ADAPTED FOR PENALTY POINTS!!!
        #print('=====================================')
        self.__icalcjac__ += 1
        p = np.array(p); n_p = len(p)        
        self.__model__.__params__.set_values(p,active_only=True) # new format, active=True
        
        #print('ENTERING RESIDUALS_JAC')
        
        # get total number of residuals and create resulting 2D matrix
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        if self.__rubber_on__ is True:
            n_tot += n_p
            
        resids = np.zeros([n_tot,n_p])
        
        # global residual index offset
        offset = 0
        
        # do "rubber" first        
        if self.__rubber_on__ is True:                  
            resids_param = self.__wpmul__*self.__model__.__wp__*np.eye(n_p) # weighted derivatives
            resids[offset:(offset+n_p)] = resids_param
            offset += n_p       
        
        #print('RESIDUALS_JAC: loop through the fit groups')
        
        # loop through the fit groups
        for grp in self.__fitgroups__:
            if not grp.__active__: continue
            #print('RESIDUALS_JAC: calculating jac on grid...')
            model_calc_jac = self.__model__.calculate_jac(grp.__inputgrid__) # calculate Jacobian, flatten the initial grid structure
            #print('RESIDUALS_JAC: ...done')
            #model_calc_jac = np.vstack(model_calc_jac) # convert array of objects to 2D array of floats; seems to be exsessive since vstack is already performed in calculate_jac
            #print('RESIDUALS_JAC: calculating func on grid...')
            model_calc = self.__model__.calculate(grp.__inputgrid__) #this step can be omitted to gain some (presumably low) speed-up
            #print('RESIDUALS_JAC: ...done')
            #print('RESIDUALS_JAC: setting calc vals and calculating residuals...')
            grp.set_calc_vals(model_calc) # set calculated values for group; can be optimized as well by "caching"
            grp.calculate_residual_derivs()
            #print('RESIDUALS_JAC: ...done')
            if self.__weighted_fit__==True:
                resids[offset:(offset+grp.length)] = model_calc_jac * np.reshape(grp.__weighted_resid_derivs__.flatten(),[grp.length,1])
            else:
                resids[offset:(offset+grp.length)] = model_calc_jac * np.reshape(grp.__unweighted_resid_derivs__.flatten(),[grp.length,1])
            offset += grp.length        
        
        # output iteration statistics
        if self.__icalcjac__%1==0:
            p_resids = p-self.__params_initial__.get_values(active_only=True)
            print('CALC JAC %5d>>>  DIST:%13.9e'%(self.__icalcjac__,np.sqrt(np.sum(p_resids**2))))

        #print('QUITTING RESIDUALS_JAC')
            
        return resids
        
    def set_weights(self,wmod):
        self.__wmod__ = np.array(wmod)
        
    #def set_penalty(self,penalty_inputgrid,penalty_output,wpenmul):
    #    self.__penalty_inputgrid__ = penalty_inputgrid
    #    self.__penalty_output__ = np.array(penalty_output)
    #    self.__wpen__ = np.ones(penalty_inputgrid.__shape__)
    #    self.__wpenmul__ = wpenmul

    def fit_minimize(self): # a family of methods from scipy.optimize.minimize
        print('USING SCIPY.OPTIMIZE.MINIMIZE: METHOD=',self.__method__)
        argv = self.__options__.copy()
        argv['bounds'] = self.__model__.__params__.get_bounds(active_only=True)
        if self.__jac__:
            #argv['jac'] = lambda p:self.residuals_jac(p)
            raise NotImplementedError # need to implement the derivatives for minimize first
        res = opt.minimize(self.sse,self.__model__.__p__,method=self.__method__,**argv)
        self.__fitgroups__.split_global_jacobian(res.jac) # split Jacobian into groups to make statistics
        return res
    
    def fit_basinhopping(self):
        print('USING SCIPY.OPTIMIZE.BASINHOPPING/ANNEAL: METHOD=',self.__method__)
        print('WARNING: PARAMETER BOUNDS ARE IGNORED')
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
        print('USING SCIPY.OPTIMIZE.LEAST_SQUARES: METHOD=',self.__method__)
        argv = self.__options__.copy()
        if self.__method__ in {'trf','dogbox'}:
            argv['bounds'] = list(zip(*self.__model__.__params__.get_bounds(active_only=True)))
        elif self.__method__ == 'lm':
            print('Warning: skipping bounds for the "lm" method')
        if self.__jac__:
            argv['jac'] = lambda p:self.residuals_jac(p)
        print('METHOD OPTIONS:',argv)
        res = opt.least_squares(self.residuals,self.__model__.__p__,
            method=self.__method__,**argv) # with user-supplied jacobian
        self.__fitgroups__.split_global_jacobian(res.jac) # split Jacobian into groups to make statistics
        self.calculate_covariance_matrix(res) # calculate covariance matrix
        return res
        
    def fit(self):        
        t = time()
        self.__params_initial__ = copy.deepcopy(self.__model__.__params__)
        icalcfun_old = self.__icalcfun__
        icalcjac_old = self.__icalcjac__
        print('BEGIN FIT')
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
        print('END FIT')
        print('%f seconds elapsed, %d func evals, %d jac evals'%\
            (time()-t,self.__icalcfun__-icalcfun_old,self.__icalcjac__-icalcjac_old))
        self.__result__ = res
        # create model final
        self.model_final = self.__model__ # create a link instead of a separate object
        
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
        
    def calculate_covariance_matrix(self,fit_result):
        """Calculate covariance CovB between parameters.
           https://stats.stackexchange.com/questions/231868/relation-between-covariance-matrix-and-jacobian-in-nonlinear-least-squares"""
        n = len(fit_result.fun)
        X = fit_result.jac
        MSE = np.sum(fit_result.fun**2)/n
        CovB = pinv(X.T @ X)*MSE
        #CovB = inv(X.T @ X)*MSE
        self.__covb__ = CovB
                
    def __repr__(self):
        if '__result__' in self.__dict__:
            return str(self.__result__)
        else:
            return '<no __result__>'
