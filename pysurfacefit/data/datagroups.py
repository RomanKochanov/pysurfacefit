import numpy as np
from numpy.linalg import pinv, inv
import jeanny3 as j

# ===============================================================================
# =========================== FIT AND PENALTY POINTS ============================
# ===============================================================================

class FitGroups:
    """
    Class to implement the groups of the fit points (penalty, outputs, etc...).
    """
    def __init__(self,grps):
        """grps - list; index - dict"""
        self.__grps__ = grps
        self.__index__ = {grp.__name__:i for i,grp in enumerate(grps)}
        
    def __getitem__(self,name):
        i = self.__index__[name]
        grp = self.__grps__[i]
        return grp
        
    def __len__(self):
        return len(self.__grps__)
        
    def __iter__(self):
        return iter(self.__grps__)     
        
    @property
    def length(self):
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        return n_tot

    def add_group(self,grp): 
        self.__index__[grp.__name__] = len(self.__grps__)
        self.__grps__.append(grp)
    
    def get_group_stats(self):
        sse_tot = 0
        stat = j.Collection()
        for grp in self.__grps__:
            if not grp.__active__: continue
            weighted_resids = grp.__weighted_resids__
            unweighted_resids = grp.__unweighted_resids__
            calc_vals = grp.__calc_vals__
            wht_sse = np.sum(weighted_resids**2)
            unwht_sse = np.sum(unweighted_resids**2)
            whts = grp.__whtmul__*grp.__weights__
            min_wht = np.min(whts); max_wht = np.max(whts)
            wht_sd = np.sqrt(wht_sse/grp.length)
            unwht_sd = np.sqrt(unwht_sse/grp.length)
            N = np.prod(weighted_resids.shape)
            min_wht_res = np.min(weighted_resids)
            max_wht_res = np.max(weighted_resids)
            min_unwht_res = np.min(unweighted_resids)
            max_unwht_res = np.max(unweighted_resids)
            #print('%30s:   N:%5d   MIN_WHT: %5.1e   MAX_WHT: %5.1e   WHT_SSE: %7.3e   UNWHT_SSE: %7.3e   WHT_SD: %7.3e   UNWHT_SD: %7.3e'%\
            #(grp.__name__,N,min_wht,max_wht,wht_sse,unwht_sse,wht_sd,unwht_sd))
            item = {
                'GROUP':grp.__name__,
                'N':N,
                'MIN_WHT':min_wht,
                'MAX_WHT':max_wht,
                'MIN_WHT_RES':min_wht_res,
                'MAX_WHT_RES':max_wht_res,
                'MIN_UNWHT_RES':min_unwht_res,
                'MAX_UNWHT_RES':max_unwht_res,
                'WHT_SSE':wht_sse,
                'UNWHT_SSE':unwht_sse,
                'WHT_SD':wht_sd,
                'UNWHT_SD':unwht_sd
            }
            stat.update(item)
        stat.order = ['GROUP','N','MIN_WHT','MAX_WHT','MIN_WHT_RES','MAX_WHT_RES',
            'MIN_UNWHT_RES','MAX_UNWHT_RES','WHT_SSE','UNWHT_SSE','WHT_SD','UNWHT_SD']
        return stat
            
#            sse_tot += wht_sse
#        print('datagroups->sse_tot>>>','%.9e'%sse_tot)
            
    def split_array(self,array,name):
        """ Split array between groups """
        offset = 0
        for grp in self.__grps__:
            if not grp.__active__: continue
            setattr(grp,name,array[offset:(offset+grp.length)])
            offset += grp.length
            
    def split_global_jacobian(self,jac):
        """
        Split "global" Jacobian matrix between the active groups.
        Needed for calculation of the fit statistics for each group separately.
        """        
        offset = 0
        for grp in self.__grps__:
            if not grp.__active__: continue
            grp.__jacobian__ = jac[offset:(offset+grp.length)]
            offset += grp.length
            
    def collect_resids(self):
        """ Collect residuals from groups for outlier statistics """
        n_tot = sum([grp.length for grp in self.__fitgroups__ if grp.__active__])
        self.__weighted_residuals__ = np.zeros(n_tot)
        self.__unweighted_residuals__ = np.zeros(n_tot)
        offset = 0
        for grp in self.__grps__:
            if not grp.__active__: continue
            self.__weighted_residuals__[offset:(offset+grp.length)] = grp.__weighted_residuals__
            self.__unweighted_residuals__[offset:(offset+grp.length)] = grp.__unweighted_residuals__
            offset += grp.length 
                   
    def calculate_outlier_statistics(self,jac):
        """Calculate useful fit statistics (cook, jackknife etc...).
           For some of the statistics Jacobian matrix is required.
           
           Used the following references:
           https://en.wikipedia.org/wiki/Leverage_(statistics)
           https://en.wikipedia.org/wiki/Studentized_residual
           https://en.wikipedia.org/wiki/DFFITS
           https://www.mathworks.com/help/stats/cooks-distance.html
           """
           
        self.collect_resids() # collect residuals
        self.load_jacobian() # load cached Jacobian

        #n = X.shape[0] # number of observables
        #m = X.shape[1] # number of free parameters
        n = self.length
        m = jac.shape[1]
           
        #X = self.__jacobian__ 
        X = jac[:n]
        
        eps = self.__weighted_resids__; eps = eps.flatten() # !!!! IF FLATTEN IS OMITTED, THEN DIMENSIONALITY MISMATCH OCCURS; ANOTHER WAY OF SOLVING THIS IS TO KEEP THE ORIGINAL DIMENSIONALITY (TODO!!!!)
        eps2 = eps**2
        sse = np.sum(eps2) # sum of squares error
        mse = sse/n # mean squared error
           
        # Calculate Hat matrix
        H = X @ pinv(X.T @ X) @ X.T # !!! using generalized inversion (pinv) since rank(X)<n_p
        
        # Leverages
        leverage = np.diag(H)
        self.__leverage__ = leverage
        
        # Studentized residuals (weighted)
        sigma2 = 1/(n-m)*sse
        T = eps/(np.sqrt(sigma2)*np.sqrt(1-leverage))
        self.__weighted_studentized_resids__ = T
        
        # DFFITS
        DFFITS = T*np.sqrt(leverage/(1-leverage))
        self.__DFFITS__ = DFFITS
        
        # Cook's distances
        cook = eps2/(m*mse)*(leverage/(1-leverage)**2)
        self.__cook__ = cook
        
        # Split statistics to groups
        self.split_array(self.__leverage__,'__leverage__')
        self.split_array(self.__weighted_studentized_resids__,'__weighted_studentized_resids__')
        self.split_array(self.__DFFITS__,'__DFFITS__')
        self.split_array(self.__cook__,'__cook__')   
        
    def calculate_covariance_matrix(self,jac):
        """Calculate covariance CovB between parameters.
           https://stats.stackexchange.com/questions/231868/relation-between-covariance-matrix-and-jacobian-in-nonlinear-least-squares"""
        n = self.length
        X = jac[:n]
        MSE = np.sum(fit_result.fun**2)/n
        CovB = pinv(X.T @ X)*MSE
        #CovB = inv(X.T @ X)*MSE
        self.__covb__ = CovB
            
    #def __repr__(self):
    #    return self.__name__

#class AbstractFitPoints(ABC):
class AbstractFitPoints:
    """
    Abstract class for Fitting points.
    Has the grids for input (N per M) and desired output (N per 1) points,
    where N is number of "observables", M is number of dimensions.
    Children must implement the residuals and their derivatives for 
    the given output values.
    """
    def __init__(self,name,input,output,whtmul=1,whts=None,active=True):
        self.__name__ = name
        self.__active__ = active
        self.__inputgrid__ = input
        self.__output__ = np.array(output)
        self.__whtmul__ = whtmul
        if whts is not None:
            self.__weights__ = np.array(whts)
        else:
            self.__weights__ = np.full(output.shape,fill_value=1.0)
        self.__calc_vals__ = None
        self.__weighted_resids__ = None
        self.__unweighted_resids__ = None
        
    def __getitem__(self,varname):
        """Get the mesh of the corresponding __inputgrid__ by the variable name"""
        return self.__inputgrid__[varname]
            
    @property
    def length(self):
        return self.__output__.size        
        
    def set_weights(self,whts):
        self.__weights__ = np.array(whts)
        
    def set_calc_vals(self,calc_vals):
        """Set last calculated values, and calculate 
           weighted and unweighted residuals"""
        self.__calc_vals__ = calc_vals
                        
    def calculate_residuals(self):
        """Array of residuals between output and calc_vals"""        
        self.__unweighted_resids__ = np.zeros(self.__output__.shape)
        # calculate unweighted resids
        for i,_ in np.ndenumerate(self.__output__):
            self.__unweighted_resids__[i] = self.resid_func(self.__output__[i],self.__calc_vals__[i])
        # calculate weighted resids
        self.__weighted_resids__ = self.__whtmul__*self.__weights__*self.__unweighted_resids__
        
    def calculate_residual_derivs(self):
        """Array of derivatives of residuals by 'calc'"""
        self.__unweighted_resid_derivs__ = np.zeros(self.__output__.shape)
        # calculate unweighted resids
        for i,_ in np.ndenumerate(self.__output__):
            self.__unweighted_resid_derivs__[i] = self.resid_deriv_func(self.__output__[i],self.__calc_vals__[i])
        # calculate weighted resids
        self.__weighted_resid_derivs__ = self.__whtmul__*self.__weights__*self.__unweighted_resid_derivs__   
        
    def stats(self,statname):
        """
        Useful shortcut to the fit statistics of all types.
        """
        if statname=='cook':
            return self.__cook__
        elif statname=='dffits':
            return self.__DFFITS__
        elif statname=='leverage':
            return self.__leverage__
        elif statname=='student':
            return self.__weighted_studentized_resids__
        else:
            raise Exception('unknown statistic %s'%statname)

    #@abstractmethod
    def resid_func(self,obs,calc):
        """Residual, scalar"""
        raise NotImplementedError

    #@abstractmethod
    def resid_deriv_func(self,obs,calc):
        """Derivative by calc, scalar"""
        raise NotImplementedError
        
    def __repr__(self):
        return self.__name__
        
class FitPoints(AbstractFitPoints): # can be sped up by re-defining resids and d_resids
    def resid_func(self,obs,calc):
        return obs-calc
    
    def resid_deriv_func(self,obs,calc):
        return -1
    
class PenaltyPointsUp(AbstractFitPoints):
    def resid_func(self,obs,calc):
        return obs-calc if obs<calc else 0
        
    def resid_deriv_func(self,obs,calc):
        return -1 if obs<calc else 0

class PenaltyPointsDown(AbstractFitPoints):
    def resid_func(self,obs,calc):
        return obs-calc if obs>calc else 0
        
    def resid_deriv_func(self,obs,calc):
        return -1 if obs>calc else 0
